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
    
    #define NUM_TISSUE  22
    
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
        return (N*sum_XY-sum_X*sum_Y) / sqrt((N*squareSum_X-sum_X*sum_X) * (N*squareSum_Y-sum_Y*sum_Y));
    }
    
    __device__ float spearman_correlation(float *X, float *Y, const int N)
    {
        float *X_rank = (float*) malloc(sizeof(float) * N);
        float *Y_rank = (float*) malloc(sizeof(float) * N);
        float coeff;
    
        rankify(X, X_rank, N);
        rankify(Y, Y_rank, N);
        coeff = correlationCoefficient(X_rank, Y_rank, N);
    
        free(X_rank);
        free(Y_rank);
        return coeff;
    }
    
    __device__ float get_sum(int start, int end, float *X_rsc_gpu, const int N)
    {
        float sum = 0;
        int rsc_start, rsc_end;
        for(int i=0; i<N; i++) {
            rsc_start = X_rsc_gpu[i * WIDTH + START];
            rsc_end = X_rsc_gpu[i * WIDTH + END];
            
            if(rsc_end < start) continue;
            else if(rsc_start > end) break;
            else {sum += X_rsc_gpu[i * WIDTH + SCORE];}
        }
        return sum;
    }
    
    __global__ void cuda_corr(float *X1, float *X2, float *out, const int N)
    {
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        if (idx >= N) return;
        
        int offset = NUM_TISSUE * idx; 
        out[idx] = correlationCoefficient(X1[offset], X2[offset], NUM_TISSUE);        
    }
    
    __global__ void cuda_sum(int *X_ref_gpu, float *X_rsc_gpu1, float *X_rsc_gpu2, float *out_buffer_gpu, int N, int M1, int M2)
    {
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        if (idx >= N) return;
    
        int ref_start = X_ref_gpu[idx * REF_WIDTH + REF_START];
        int ref_end = X_ref_gpu[idx * REF_WIDTH + REF_END];
        
        out_buffer_gpu[idx * OUT_WIDTH + SUM1] = get_sum(ref_start, ref_end, X_rsc_gpu1, M1);
        out_buffer_gpu[idx * OUT_WIDTH + SUM2] = get_sum(ref_start, ref_end, X_rsc_gpu2, M2);
    }""")


class Correlation:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.table_names = {}

    def to_server(self):
        which = 'newton'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='04:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def correlation_fan_rna_cpu(self):
        # tissue list
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'tissues_fantom_rna.xlsx')
        df_tis = pd.read_excel(fpath, sheet_name='Sheet2')
        tissues = df_tis['RNA-seq'].iloc[:22]

        # reference GENCODE
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation_spt.db')
        ref_con = sqlite3.connect(ref_path)

        # Fantom5
        fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'fantom_cage_by_tissue_spt.db')
        con_fan = sqlite3.connect(fan_path)

        # RNA-seq
        rna_path = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_spt.db')
        con_rna = sqlite3.connect(rna_path)

        # out
        M = len(tissues)
        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna.db')
        con_out = sqlite3.connect(out_path)

        columns = ['chromosome', 'start', 'end', 'strand', 'gene_name', 'transcript_name', 'transcript_type', 'corr']
        tlist = Database.load_tableList(ref_con)
        dfs = []
        for tname in tlist:
            source, version, chromosome, strand = tname.split('.')
            df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)
            N = df_ref.shape[0]
            print(chromosome, strand)

            df_res = pd.DataFrame(index=df_ref.index, columns=columns)
            for idx in df_ref.index:
                if idx % 100 == 0 or idx + 1 == df_ref.shape[0]:
                    print('{:,d} / {:,d}'.format(idx + 1, df_ref.shape[0]))
                if strand == '+':
                    tss = df_ref['start']
                else:
                    tss = df_ref['end']

                df_ref['end'] = tss + 500
                df_ref['start'] = tss - 500

                start = df_ref.loc[idx, 'start']
                end = df_ref.loc[idx, 'end']
                gene_name = df_ref.loc[idx, 'gene_name']
                transcript_name = df_ref.loc[idx, 'transcript_name']
                transcript_type = df_ref.loc[idx, 'transcript_type']

                buffer = np.zeros((M, 2))
                for i, tissue in enumerate(tissues):
                    tname_rsc = '_'.join([tissue, chromosome, strand])
                    df_fan = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND strand='{}' AND NOT "
                                               "start>{end} AND NOT end<{start}".format(tname_rsc, chromosome, strand,
                                                                                        start=start, end=end), con_fan)
                    df_rna = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND strand='{}' AND NOT "
                                               "start>{end} AND NOT end<{start}".format(tname_rsc, chromosome, strand,
                                                                                        start=start, end=end), con_rna)
                    buffer[i, 0] = df_fan['score'].sum()
                    buffer[i, 1] = df_rna['FPKM'].sum()

                corr = np.corrcoef(buffer[:, 0], buffer[:, 1])
                df_res.loc[idx, :] = [chromosome, start, end, strand, gene_name, transcript_name, transcript_type, corr[0, 1]]
            dfs.append(df_res)
        pd.concat(dfs).to_sql(tname, con_out, if_exists='replace', index=None)

    def set_gpu_buffer(self, df_ref, df_rsc):
        N = np.int32(df_ref.shape[0])

        X_ref = df_ref.values.flatten().astype(np.int32)
        X_ref_gpu = cuda.mem_alloc(X_ref.nbytes)
        cuda.memcpy_htod(X_ref_gpu, X_ref)

        M = []
        X_rsc_gpu = []
        for df in df_rsc:
            M.append(np.int32(df.shape[0]))
            rsc = df.values.flatten().astype(np.float32)
            rsc_gpu = cuda.mem_alloc(rsc.nbytes)
            cuda.memcpy_htod(rsc_gpu, rsc)
            X_rsc_gpu.append(rsc_gpu)

        return X_ref_gpu, X_rsc_gpu, N, M

    def cuda_sum(self, X_ref_gpu, X_rsc_gpu, N, M):
        THREADS_PER_BLOCK = 1 << 10
        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)

        out_buffer = np.zeros((N, 2)).flatten().astype(np.float32)
        out_buffer_gpu = cuda.mem_alloc(out_buffer.nbytes)

        func = mod.get_function("cuda_sum")
        func(X_ref_gpu, X_rsc_gpu[0], X_rsc_gpu[1], out_buffer_gpu, N, M[0], M[1], block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))

        cuda.memcpy_dtoh(out_buffer, out_buffer_gpu)
        out_buffer = out_buffer.reshape((N, 2))
        return out_buffer

    def cuda_corr(self, X, Y, N):
        THREADS_PER_BLOCK = 1 << 10
        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)

        X = X.values.flatten().astype(np.float32)
        Y = Y.values.flatten().astype(np.float32)
        X_gpu = cuda.mem_alloc(X.nbytes)
        Y_gpu = cuda.mem_alloc(Y.nbytes)

        cuda.memcpy_htod(X_gpu, X)
        cuda.memcpy_htod(Y_gpu, Y)

        out = np.zeros(N).astype(np.float32)
        out_gpu = cuda.mem_alloc(out.nbytes)

        func = mod.get_function("cuda_corr")
        func(X, Y, out_gpu, N, block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out, out_gpu)
        return out

    def correlation_fan_rna(self):
        # tissue list
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'tissues_fantom_rna.xlsx')
        df_tis = pd.read_excel(fpath, sheet_name='Sheet2')
        tissues = df_tis['RNA-seq'].iloc[:22]

        # reference GENCODE
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation_spt.db')
        ref_con = sqlite3.connect(ref_path)

        # Fantom5
        fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'fantom_cage_by_tissue_spt.db')
        con_fan = sqlite3.connect(fan_path)

        # RNA-seq
        rna_path = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_spt.db')
        con_rna = sqlite3.connect(rna_path)

        # output
        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna.db')
        con_out = sqlite3.connect(out_path)

        M_ = len(tissues)

        tlist = Database.load_tableList(ref_con)
        for tname in tlist:
            source, version, chromosome, strand = tname.split('.')
            df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)
            if df_ref.empty:
                continue

            if strand == '+':
                tss = deepcopy(df_ref['start'])
            else:
                tss = deepcopy(df_ref['end'])

            df_ref['start'] = tss - 500
            df_ref['end'] = tss + 500

            N_ = df_ref.shape[0]
            print(chromosome, strand)

            buffer = np.zeros((N_, 2, M_))
            for i, tissue in enumerate(tissues):
                tname_rsc = '_'.join([tissue, chromosome, strand])
                df_fan = pd.read_sql_query("SELECT start, end, score FROM '{}' WHERE chromosome='{}' AND strand='{}'"
                                           "".format(tname_rsc, chromosome, strand), con_fan).sort_values(by=['start'])
                df_rna = pd.read_sql_query("SELECT start, end, FPKM FROM '{}' WHERE chromosome='{}' AND strand='{}'"
                                           "".format(tname_rsc, chromosome, strand), con_rna).sort_values(by=['start'])
                if df_fan.empty:
                    buffer[:, 0, i] = 0
                elif df_rna.empty:
                    buffer[:, 1, i] = 0
                else:
                    X_ref_gpu, X_rsc_gpu, N, M = self.set_gpu_buffer(df_ref[['start', 'end']], [df_fan, df_rna])
                    buffer[:, :, i] = self.cuda_sum(X_ref_gpu, X_rsc_gpu, N, M)
            
            df_ref['corr'] = self.cuda_corr(buffer[:, 0, :], buffer[:, 0, :], N_)
            df_ref.to_sql(tname, con_out, index=None, if_exists='replace')


if __name__ == '__main__':
    cor = Correlation()
    if cor.hostname == 'mingyu-Precision-Tower-7810':
        cor.to_server()
    else:
        cor.correlation_fan_rna()
        # cor.run()
