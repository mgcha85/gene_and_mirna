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

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
from Database import Database

mod = SourceModule("""
#include <stdio.h>
enum{TABLE_NUM, FPKM_OUT, TPM_OUT, OUT_WIDTH};
enum{START, END, FPKM, TPM, WIDTH};
enum{REF_START, REF_END, REF_WIDTH};

__device__ float get_rpkm(int *res_buffer, const int ref_start, const int ref_end, const int N, const int idx)
{
    int cnt=0;
    int res_start, res_end;
    int length = ref_end - ref_start;
    float rpkm = 0.0;

    for(int i=0; i<N; i++) {
        res_start = res_buffer[i * WIDTH + START];      
        res_end = res_buffer[i * WIDTH + END];      
        if(ref_start > res_end || ref_end < res_start)  continue;
        else                                            cnt++;
    }
    rpkm = (cnt / 1e6) / length;
    return rpkm;
}

__device__ void get_reads_count(float *res_buffer, float* out_buffer_gpu, const int ref_start, const int ref_end, const int N, const int idx)
{
    int cnt=0;
    int res_start, res_end;
    float fpkm=0.0, tpm=0.0;
    
    for(int i=0; i<N; i++) {
        res_start = res_buffer[i * WIDTH + START];      
        res_end = res_buffer[i * WIDTH + END];  
        if(ref_start > res_end || ref_end < res_start)  
            continue;
        else {
            fpkm += res_buffer[i * WIDTH + FPKM];
            tpm += res_buffer[i * WIDTH + TPM];
        }
    }
    out_buffer_gpu[idx * OUT_WIDTH + FPKM_OUT] = fpkm;
    out_buffer_gpu[idx * OUT_WIDTH + TPM_OUT] = tpm;
}

__device__ int get_table_num(int *ref_data_lengths_cum_gpu, const int idx, const int N)
{
    for(int i=0; i<N-1; i++) {  
        if(ref_data_lengths_cum_gpu[i] <= idx && ref_data_lengths_cum_gpu[i+1] > idx) return i;
    }
    return N;
}

__global__ void cuda_scanner(int *ref_buffer_gpu, int *ref_data_lengths_cum_gpu, float *res_buffer_gpu, int *res_data_lengths_cum_gpu, float *out_buffer_gpu, const int N, const int M)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N) return;

    int tb_num, sidx, eidx;
    int ref_start, ref_end;

    ref_start = ref_buffer_gpu[idx * REF_WIDTH + REF_START];
    ref_end = ref_buffer_gpu[idx * REF_WIDTH + REF_END];

    tb_num = get_table_num(ref_data_lengths_cum_gpu, idx, M);

    sidx = res_data_lengths_cum_gpu[tb_num];
    eidx = res_data_lengths_cum_gpu[tb_num + 1];

    out_buffer_gpu[idx * OUT_WIDTH + TABLE_NUM] = (float) tb_num;
    get_reads_count(&res_buffer_gpu[sidx * WIDTH], out_buffer_gpu, ref_start, ref_end, eidx - sidx, idx);

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
        server = Server()
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server.job_script(fname)

        server_root = os.path.join(server.server, 'source/gene_and_mirna')
        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def split_table_for_res(self, df):
        df_sorted = []
        dfs = []
        data_length = []
        df_chr = df.groupby('chromosome')

        for tname, tnum in self.table_names.items():
            chr, str = tname.split(';')
            if len(chr) > 5 or chr not in df_chr.groups:
                continue
            df_sub = df_chr.get_group(chr)
            # df_sub_sub = df_sub[df_sub['strand'] == str]

            dfs.append(df_sub[['start', 'end', 'FPKM', 'TPM']].values.flatten().astype(np.float32))
            data_length.append(df_sub.shape[0])
            df_sorted.append(df_sub)
        return dfs, data_length, pd.concat(df_sorted).reset_index(drop=True)

    def split_table_for_ref(self, df):
        RANGE = 500

        df_sorted = []
        dfs = []
        data_length = []
        df_chr = df.groupby('chromosome')
        chromosomes = sorted(df_chr.groups)

        cnt = 0
        for chr in chromosomes:
            if len(chr) > 5:
                continue
            df_sub = df_chr.get_group(chr)
            df_sub_str = df_sub.groupby('strand')
            strands = sorted(df_sub_str.groups)
            for str in strands:
                tname = '{};{}'.format(chr, str)
                self.table_names[tname] = cnt
                cnt += 1
                df_sub_sub = df_sub_str.get_group(str)
                if str == '+':
                    tss = deepcopy(df_sub_sub['start'])
                else:
                    tss = deepcopy(df_sub_sub['end'])

                df_sub_sub.loc[:, 'start'] = tss - RANGE
                df_sub_sub.loc[:, 'end'] = tss + RANGE
                dfs.append(df_sub_sub[['start', 'end']].values.flatten().astype(np.int32))
                data_length.append(df_sub_sub.shape[0])
                df_sorted.append(df_sub_sub)
        return dfs, data_length, pd.concat(df_sorted).reset_index(drop=True)

    def set_data(self, dfs, data_length):
        dfs = np.concatenate(dfs)
        data_length = np.array(data_length)
        data_length_cum = np.zeros(data_length.shape[0] + 1)
        data_length_cum[1:] = data_length.cumsum()
        return dfs, data_length_cum.astype(np.int32)

    def load_reference(self, tissue):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.final_confirmed_tss.db')
        con = sqlite3.connect(fpath)
        if Database.checkTableExists(con, tissue):
            return pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con)

    def to_gpu(self, dfs_ref, data_length_cum_ref, dfs_res, data_length_cum_src, N):
        THREADS_PER_BLOCK = 1 << 10
        OUT_WDITH = 3

        ref_buffer_gpu = cuda.mem_alloc(dfs_ref.nbytes)
        cuda.memcpy_htod(ref_buffer_gpu, dfs_ref)

        res_buffer_gpu = cuda.mem_alloc(dfs_res.nbytes)
        cuda.memcpy_htod(res_buffer_gpu, dfs_res)

        ref_data_lengths_cum_gpu = cuda.mem_alloc(data_length_cum_ref.nbytes)
        cuda.memcpy_htod(ref_data_lengths_cum_gpu, data_length_cum_ref)

        res_data_lengths_cum_gpu = cuda.mem_alloc(data_length_cum_src.nbytes)
        cuda.memcpy_htod(res_data_lengths_cum_gpu, data_length_cum_src)

        M = np.int32(data_length_cum_ref.shape[0])

        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)
        out_buffer = np.zeros((N, OUT_WDITH)).flatten().astype(np.float32)
        out_buffer_gpu = cuda.mem_alloc(out_buffer.nbytes)

        func = mod.get_function("cuda_scanner")
        func(ref_buffer_gpu, ref_data_lengths_cum_gpu, res_buffer_gpu, res_data_lengths_cum_gpu, out_buffer_gpu, N, M,
             block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out_buffer, out_buffer_gpu)
        out_buffer = out_buffer.reshape((-1, OUT_WDITH))
        return pd.DataFrame(data=out_buffer, columns=['TABLE_NUM', 'FPKM', 'TPM'])

    def get_tpm(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue.db')
        con_fan = sqlite3.connect(fpath)

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con_rna = sqlite3.connect(fpath)

        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'correlation_fantom.db')
        con_out = sqlite3.connect(out_path)

        tlist_fan = Database.load_tableList(con_fan)
        tlist_rna = Database.load_tableList(con_rna)

        for tname_fan in tlist_fan:
            print(tname_fan)
            df_ref = self.load_reference(tname_fan)
            if df_ref is None:
                continue

            dfs_ref, data_length, df_sorted = self.split_table_for_ref(df_ref)
            dfs_ref, data_length_cum_ref = self.set_data(dfs_ref, data_length)

            df_fan = pd.read_sql_query("SELECT * FROM {}".format(tname_fan), con_fan)
            dfs_fan, data_length_fan, df_sorted_fan = self.split_table_for_res(df_fan)
            dfs_fan, data_length_cum_src = self.set_data(dfs_fan, data_length_fan)
            df_out = self.to_gpu(dfs_ref, data_length_cum_ref, dfs_fan, data_length_cum_src, np.int32(df_sorted.shape[0]))
            df_out = pd.concat([df_sorted, df_out], axis=1)

            df_out.to_sql(tname_fan, con_out, if_exists='replace', index=None)

            tnames_rna = [x for x in tlist_rna if tname_fan in x]
            for tname_rna in tnames_rna:
                df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tname_rna), con_rna)
                dfs_rna, data_length_rna, df_sorted_rna = self.split_table_for_res(df_rna)
                dfs_rna, data_length_cum_src = self.set_data(dfs_rna, data_length_rna)

                df_out = self.to_gpu(dfs_ref, data_length_cum_ref, dfs_rna, data_length_cum_src, np.int32(df_ref.shape[0]))
                df_out = pd.concat([df_sorted, df_out], axis=1)
                df_out.to_sql(tname_rna, con_out, if_exists='replace', index=None)

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'correlation_rna_seq.db')
        con_rna = sqlite3.connect(fpath)
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'correlation_fantom.db')
        con_fan = sqlite3.connect(fpath)

        tlist_fan = Database.load_tableList(con_fan)
        tlist_rna = Database.load_tableList(con_rna)

        report = []
        for tname in tlist_fan:
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_fan)

            tlist = [x for x in tlist_rna if tname in x]
            for tname_rna in tlist:
                df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tname_rna), con_rna)
                fpkm = scipy.stats.spearmanr(df_fan['FPKM'], df_rna['FPKM'])
                tpm = scipy.stats.spearmanr(df_fan['TPM'], df_rna['TPM'])
                report.append([tname, tname_rna, fpkm, tpm])


if __name__ == '__main__':
    cor = Correlation()
    if cor.hostname == 'mingyu-Precision-Tower-7810':
        cor.to_server()
    else:
        cor.get_tpm()