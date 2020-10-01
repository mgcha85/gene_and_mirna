import pandas as pd
from Database import Database
import numpy as np
from itertools import product
import socket

import os
import sqlite3

import pyopencl as cl

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

mf = cl.mem_flags

prg = cl.Program(ctx, """
void rankify(__global double *X, __global double *Rank_X, const int N) {
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

double correlationCoefficient(__global const double *X, __global const double *Y, const int N)
{
    double sum_X=0, sum_Y=0, sum_XY=0;
    double squareSum_X=0, squareSum_Y=0;

    for(int i=0; i<N; i++)
    {
        sum_X = sum_X + X[i];
        sum_Y = sum_Y + Y[i];
        sum_XY = sum_XY + X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X = squareSum_X + X[i] * X[i];
        squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }
    double den = sqrt((N*squareSum_X - sum_X*sum_X) * (N*squareSum_Y - sum_Y*sum_Y));
    double nom = (N*sum_XY - sum_X*sum_Y);
    return nom / den;
}

double spearman_correlation(__global double *X, __global double *Y, __global double *X_rank, 
                            __global double *Y_rank, const int N)
{
    rankify(X, X_rank, N);
    rankify(Y, Y_rank, N);
    return correlationCoefficient(X_rank, Y_rank, N);
}

__kernel void cl_spearman(__global double *X, __global double *Y, __global double *X_rank, 
                            __global double *Y_rank, __global double *out, __constant int *nums)
{
    int idx = get_global_id(0);
    
    int N = nums[0];
    int M = nums[1];
    int WIDTH = nums[2];
    int prod = nums[3];
    
    if(idx == 0)
        printf("N: %d, M: %d, WIDTH: %d, prod: %d\\n", N, M, WIDTH, prod);
    
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
    out[idx] = spearman_correlation(&X[nidx], &Y[midx], &X_rank[nidx], &Y_rank[midx], WIDTH);
}

__kernel void cl_pearson(__global double *X, __global double *Y, __global double *out, __constant double *nums)
{
    int idx = get_global_id(0);
    int N = nums[0];
    int M = nums[1];
    int WIDTH = nums[2];
    int prod = nums[3];

    int NM, nidx, midx;
    if(prod == 1)   NM = N * M;
    else            NM = N;

    if (idx >= NM) return;

    if(prod == 1) {
        nidx = (idx / M) * WIDTH;
        midx = (idx % M) * WIDTH;
    }
    else {
        nidx = idx * WIDTH;
        midx = idx * WIDTH;
    }
    out[idx] = correlationCoefficient(&X[nidx], &Y[midx], WIDTH);
}
""").build()



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

        N = np.int32(X.shape[0])
        M = np.int32(Y.shape[0])
        WIDTH = np.int32(X.shape[1])
        p = np.int32(prod)
        if prod:
            NM = N * M
        else:
            NM = N

        gene = X.values.flatten().astype(np.float64)
        gene_gpu = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gene)

        gene_rank = np.zeros_like(gene).astype(np.float64)
        gene_rank_gpu = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gene_rank)

        mir = Y.values.flatten().astype(np.float64)
        mir_gpu = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=mir)

        mir_rank = np.zeros_like(mir).astype(np.float64)
        mir_rank_gpu = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=mir_rank)

        out = np.zeros(NM).astype(np.float64)
        out_gpu = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=out)

        nums = np.array([N, M, WIDTH, p]).astype(np.int)
        nums_gpu = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=nums)

        prg.cl_spearman(queue, (NM, ), None, gene_gpu, mir_gpu, gene_rank_gpu, mir_rank_gpu, out_gpu, nums_gpu)
        cl.enqueue_copy(queue, out, out_gpu)
        df_res = pd.DataFrame(out.reshape((N, M)), index=X.index, columns=Y.index)
        return df_res

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

    def test_cpu(self, hbw):
        from scipy.stats import spearmanr
        from joblib import Parallel, delayed
        import multiprocessing

        def processInput(i, gidx):
            if (i + 1) % 100 == 0 or (i+1) == df_gene.shape[0]:
                print('{} / {}'.format(i+1, df_gene.shape[0]))
            df_out = pd.Series(index=df_mir.index)
            for midx in df_mir.index:
                df_out[midx] = spearmanr(df_gene.loc[gidx, '10964C':], df_mir.loc[midx, '10964C':])[0]
            return df_out

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
        num_cores = multiprocessing.cpu_count()

        dfs = Parallel(n_jobs=num_cores)(delayed(processInput)(i, gidx) for i, gidx in enumerate(df_gene.index))
        df_out = pd.concat(dfs, axis=1).T
        df_out.index = df_gene.index
        df_out.to_sql('corr_cpu', con_out, if_exists='replace')


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
        N = np.int32(X.shape[0])
        M = np.int32(Y.shape[0])
        WIDTH = np.int32(X.shape[1])
        p = np.int32(prod)
        if prod:
            NM = N * M
        else:
            NM = N

        gene = X.values.flatten().astype(np.float32)
        gene_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gene)

        mir = Y.values.flatten().astype(np.float32)
        mir_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=mir)

        mir_rank = np.zeros_like(Y).flatten().astype(np.float32)
        mir_rank_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=mir_rank)

        out = np.zeros(NM).astype(np.float32)
        out_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=out)

        prg.cl_pearson(queue, (N, M), None, gene_g, mir_g, out_g)
        cl.enqueue_copy(queue, out, out_g)
        df_res = pd.DataFrame(out.reshape((N, M)), index=X.index, columns=Y.index)
        return df_res

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


if __name__ == '__main__':
    hostname = socket.gethostname()
    if hostname == 'mingyu-Precision-Tower-7810':
        root = '/home/mingyu/Bioinformatics'
    elif hostname == 'DESKTOP-DLOOJR6' or hostname == 'DESKTOP-1NLOLK4':
        root = 'D:/Bioinformatics'
    elif hostname == 'mingyu-Inspiron-7559':
        root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
    else:
        root = '/lustre/fs0/home/mcha/Bioinformatics'

    spr = Spearman(root)
    spr.test_cpu(100)
