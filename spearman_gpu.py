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

    enum{REF_START, REF_END, REF_WIDTH};
    enum{START, END, SCORE, WIDTH};
    enum{SUM1, SUM2, OUT_WIDTH};

    #define NUM_CELLLINE  240

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
        return (N*sum_XY - sum_X*sum_Y) / sqrt((N*squareSum_X - sum_X*sum_X) * (N*squareSum_Y - sum_Y*sum_Y));
    }

    __device__ float spearman_correlation(float *X, float *Y, float *X_rank, float *Y_rank, const int N)
    {
        float coeff;

        rankify(X, X_rank, N);
        rankify(Y, Y_rank, N);
        coeff = correlationCoefficient(X_rank, Y_rank, N);

        free(X_rank);
        free(Y_rank);
        return coeff;
    }""")


class Spearman:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.table_names = {}

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

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def run(self, X, Y, df_out):
        # X, Y: pandas DataFrame

        comb = np.array(list((product(range(X.shape[0]), range(Y.shape[0])))))

        # gpu
        THREADS_PER_BLOCK = 1 << 10
        N = X.shape[0]
        gene = X.values.flatten().astype(np.int32)
        gene_gpu = cuda.mem_alloc(gene.nbytes)
        cuda.memcpy_htod(gene_gpu, gene)

        M = Y.shape[0]
        mir = Y.values.flatten().astype(np.int32)
        mir_gpu = cuda.mem_alloc(mir.nbytes)
        cuda.memcpy_htod(mir_gpu, mir)

        idx = comb.flatten().astype(np.int32)
        idx_gpu = cuda.mem_alloc(idx.nbytes)
        cuda.memcpy_htod(idx_gpu, idx)

        out = np.zeros(df_out.shape)
        out_gpu = cuda.mem_alloc(out.nbytes)
        cuda.memcpy_htod(out_gpu, out)

        gridN = int((N * M + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)
        func = mod.get_function("cuda_corr")
        func(gene_gpu, mir_gpu, idx_gpu, out_gpu, N, M, block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out, out_gpu)
        df_out.values = out
        return df_out

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
        df_out = pd.DataFrame(index=df_gene.index, columns=df_mir.index)
        df_out = self.run(df_gene, df_mir, df_out)
        df_out.to_sql('corr', con_out, if_exists='replace')


if __name__ == '__main__':
    cor = Spearman()
    cor.test(100)