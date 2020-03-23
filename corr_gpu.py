import socket
import os
import sqlite3
import pandas as pd
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import scipy.stats
from Server import Server
import sys
from Database import Database
from copy import deepcopy
from joblib import Parallel, delayed
import multiprocessing
from scipy.stats import spearmanr
from itertools import product

hostname = socket.gethostname()
if hostname != 'mingyu-Precision-Tower-7810':
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule

    mod = SourceModule("""
    #include <stdio.h>
    #include <math.h>

    enum{NIDX, MIDX, NUM_IDX};

    __device__ void rankify(float *X, float *Rank_X, const int N) {
        for(int i=0; i<N; i++)
        {
            int r=1, s=1;

            for(int j=0; j<i; j++) {
                if(X[j] < X[i]) r++;
                if(X[j] == X[i]) s++;
            }
            for(int j=i+1; j<N; j++) {
                if (X[j] < X[i]) r++;
                if (X[j] == X[i]) s++;
            }
            Rank_X[i] = r + (s-1) * 0.5;
        }
    }

    __device__ float correlationCoefficient(float *X, float *Y, const int N)
    {
        float sum_X=0, sum_Y=0, sum_XY=0;
        float squareSum_X=0, squareSum_Y=0;

        for(int i=0; i<N; i++)
        {
            sum_X = sum_X + X[i];
            sum_Y = sum_Y + Y[i];
            sum_XY = sum_XY + X[i] * Y[i];

            // sum of square of array elements.
            squareSum_X = squareSum_X + X[i] * X[i];
            squareSum_Y = squareSum_Y + Y[i] * Y[i];
        }
        float den = sqrt((N*squareSum_X - sum_X*sum_X) * (N*squareSum_Y - sum_Y*sum_Y));
        float nom = (N*sum_XY - sum_X*sum_Y);
        //printf("den: %f, nom: %d\\n", den, nom);
        return nom / den;
    }

    __device__ float spearman_correlation(float *X, float *Y, float *X_rank, float *Y_rank, const int N)
    {
        rankify(X, X_rank, N);
        rankify(Y, Y_rank, N);
        return correlationCoefficient(X_rank, Y_rank, N);
    }

    __device__ void print_vector(float *X, const int N)
    {
        for(int i=0; i<N; i++)
            printf("%f\\n", X[i]);
    }

    __global__ void cuda_spearman(float *X, float *Y, float *X_rank, float *Y_rank, float *out, const int N, const int M, const int WIDTH, const int prod)
    {
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        int NM, nidx, midx;
        if(prod == 1)   NM = N * M;
        else            NM = N;
        
        if (idx >= NM) return;
        
        if(prod == 1) {
            nidx = (int)(idx / M) * WIDTH;
            midx = (idx % M) * WIDTH;
        }
        else {
            nidx = idx * WIDTH;
            midx = idx * WIDTH;
        }
        //printf("nidx: %d, midx: %d\\n", (int)(idx / M), (idx % M));
        out[idx] = spearman_correlation(&X[nidx], &Y[midx], &X_rank[nidx], &Y_rank[midx], WIDTH);
    }

    __global__ void cuda_pearson(float *X, float *Y, float *out, const int N, const int M, const int WIDTH, const int prod)
    {
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        int NM, nidx, midx;
        if(prod == 1)   NM = N * M;
        else            NM = N;

        if (idx >= NM) return;

        if(prod == 1) {
            nidx = (int)(idx / M) * WIDTH;
            midx = (idx % M) * WIDTH;
        }
        else {
            nidx = idx * WIDTH;
            midx = idx * WIDTH;
        }
        out[idx] = correlationCoefficient(&X[nidx], &Y[midx], WIDTH);
    }""")


class Spearman:
    def __init__(self, root):
        self.root = root

    def to_server(self):
        from Server import Server
        import sys

        which = 'newton'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='08:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command(
            "cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def run(self, X, Y, prod=True):
        # X, Y: pandas DataFrame
        # gpu
        # 2 which is not in -1 to 1 means nan

        THREADS_PER_BLOCK = 1 << 10

        # print(X.shape)
        # print(Y.shape)
        # X = X[X.sum(axis=1) > 0]
        # Y = Y[Y.sum(axis=1) > 0]
        # print(X.shape)
        # print(Y.shape)

        N = np.int32(X.shape[0])
        M = np.int32(Y.shape[0])
        WIDTH = np.int32(X.shape[1])
        p = np.int32(prod)
        if prod:
            NM = N * M
        else:
            NM = N

        gene = X.values.flatten().astype(np.float32)
        gene_gpu = cuda.mem_alloc(gene.nbytes)
        cuda.memcpy_htod(gene_gpu, gene)

        gene_rank_gpu = cuda.mem_alloc(gene.nbytes)
        cuda.memcpy_htod(gene_rank_gpu, gene)

        mir = Y.values.flatten().astype(np.float32)
        mir_gpu = cuda.mem_alloc(mir.nbytes)
        cuda.memcpy_htod(mir_gpu, mir)

        mir_rank_gpu = cuda.mem_alloc(mir.nbytes)
        cuda.memcpy_htod(mir_rank_gpu, mir)

        out = np.zeros(NM).astype(np.float32)
        out_gpu = cuda.mem_alloc(out.nbytes)
        cuda.memcpy_htod(out_gpu, out)

        gridN = int((NM + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)
        func = mod.get_function("cuda_spearman")
        func(gene_gpu, mir_gpu, gene_rank_gpu, mir_rank_gpu, out_gpu, N, M, WIDTH, p, block=(THREADS_PER_BLOCK, 1, 1),
             grid=(gridN, 1))

        cuda.memcpy_dtoh(out, out_gpu)
        if prod:
            out = out.reshape((N, M))
        return out

    def test(self, hbw):
        fpath_mir = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_mir_{}.db'.format(hbw))
        con_mir = sqlite3.connect(fpath_mir)

        fpath_gene = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_gene_{}.db'.format(hbw))
        con_gene = sqlite3.connect(fpath_gene)

        fpath_out = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_corr_{}.db'.format(hbw))
        con_out = sqlite3.connect(fpath_out)

        def merge(con):
            tlist = Database.load_tableList(con)
            dfs = []
            for tname in tlist:
                dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
            return pd.concat(dfs).set_index(dfs[-1].columns[0])

        df_mir = merge(con_mir)
        df_gene = merge(con_gene)
        df_out = pd.DataFrame(data=self.run(df_gene.loc[:, '10964C':], df_mir.loc[:, '10964C':]), index=df_gene.index,
                              columns=df_mir.index)
        df_out.to_sql('corr', con_out, if_exists='replace')


class Pearson:
    def __init__(self, root):
        self.root = root

    def to_server(self):
        from Server import Server
        import sys

        which = 'newton'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='08:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command(
            "cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def run(self, X, Y, prod=True):
        # X, Y: pandas DataFrame
        if prod:
            comb = np.array(list((product(range(X.shape[0]), range(Y.shape[0])))))
        else:
            comb = np.array(0)

        # gpu
        THREADS_PER_BLOCK = 1 << 10
        N = np.int32(X.shape[0])
        M = np.int32(Y.shape[0])
        WIDTH = np.int32(X.shape[1])
        p = np.int32(prod)
        if prod:
            NM = N * M
        else:
            NM = N

        gene = X.values.flatten().astype(np.float32)
        gene_gpu = cuda.mem_alloc(gene.nbytes)
        cuda.memcpy_htod(gene_gpu, gene)

        gene_rank_gpu = cuda.mem_alloc(gene.nbytes)
        cuda.memcpy_htod(gene_rank_gpu, gene)

        mir = Y.values.flatten().astype(np.float32)
        mir_gpu = cuda.mem_alloc(mir.nbytes)
        cuda.memcpy_htod(mir_gpu, mir)

        mir_rank_gpu = cuda.mem_alloc(mir.nbytes)
        cuda.memcpy_htod(mir_rank_gpu, mir)

        idx = comb.flatten().astype(np.int32)
        idx_gpu = cuda.mem_alloc(idx.nbytes)
        cuda.memcpy_htod(idx_gpu, idx)

        out = np.zeros(NM).astype(np.float32)
        out_gpu = cuda.mem_alloc(out.nbytes)
        cuda.memcpy_htod(out_gpu, out)

        gridN = int((NM + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)
        func = mod.get_function("cuda_pearson")
        func(gene_gpu, mir_gpu, gene_rank_gpu, mir_rank_gpu, idx_gpu, out_gpu, N, M, WIDTH, p, block=(THREADS_PER_BLOCK, 1, 1),
             grid=(gridN, 1))

        cuda.memcpy_dtoh(out, out_gpu)
        if prod:
            out = out.reshape((N, M))
        return out

    def test(self, hbw):
        fpath_mir = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_mir_{}.db'.format(hbw))
        con_mir = sqlite3.connect(fpath_mir)

        fpath_gene = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_gene_{}.db'.format(hbw))
        con_gene = sqlite3.connect(fpath_gene)

        fpath_out = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_corr_{}.db'.format(hbw))
        con_out = sqlite3.connect(fpath_out)

        def merge(con):
            tlist = Database.load_tableList(con)
            dfs = []
            for tname in tlist:
                dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
            return pd.concat(dfs).set_index(dfs[-1].columns[0])

        df_mir = merge(con_mir)
        df_gene = merge(con_gene)
        df_out = pd.DataFrame(data=self.run(df_gene.loc[:, '10964C':], df_mir.loc[:, '10964C':]), index=df_gene.index,
                              columns=df_mir.index)
        df_out.to_sql('corr', con_out, if_exists='replace')
