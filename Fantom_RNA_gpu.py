import os
import sqlite3
import pandas as pd
import numpy as np
import socket
from Database import Database
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
from Database import Database

mod = SourceModule("""
#include <stdio.h>
enum{TABLE_NUM, RPKM, OUT_WIDTH};
enum{START, END, WIDTH};

__device__ int get_rpkm(int *res_buffer, const int ref_start, const int ref_end, const int N)
{
    int cnt=0, rpkm=0;
    int res_start, res_end;
    int length = ref_end - ref_start;
    
    for(int i=0; i<N; i++) {
        res_start = res_buffer[i * WIDTH + START];      
        res_end = res_buffer[i * WIDTH + END];      
        if(ref_start > res_end || ref_end < res_start)  continue;
        else                                            cnt++;
    }
    rpkm = (cnt / 1e6) / length;
    return rpkm;
}

__device__ int get_table_num(int *ref_data_lengths_cum_gpu, const int idx, const int N)
{
    for(int i=0; i<N; i++) {  
        if(ref_data_lengths_cum_gpu[i] <= idx && ref_data_lengths_cum_gpu[i+1] > idx) return i;
    }
    return N;
}

__global__ void cuda_scanner(int *ref_buffer_gpu, int *ref_data_lengths_cum_gpu, int *res_buffer_gpu, int *res_data_lengths_cum_gpu, int *out_buffer_gpu, const int N, const int M)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N) return;

    int tb_num, sidx, eidx;
    int ref_start, ref_end;
    
    ref_start = ref_buffer_gpu[idx * WIDTH + START];
    ref_end = ref_buffer_gpu[idx * WIDTH + END];
    
    tb_num = get_table_num(ref_data_lengths_cum_gpu, idx, M);
    
    sidx = res_data_lengths_cum_gpu[tb_num];
    eidx = res_data_lengths_cum_gpu[tb_num + 1];

    out_buffer_gpu[idx * OUT_WIDTH + TABLE_NUM] = tb_num;
    out_buffer_gpu[idx * OUT_WIDTH + RPKM] = get_rpkm(&res_buffer_gpu[sidx * WIDTH], ref_start, ref_end, idx, eidx - sidx);

}""")


class Fantom_RNA:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.cells = []

    def tpm(self, df_gene, df_rna):
        length = abs(df_gene['start'] - df_gene['end'])
        N = df_rna.shape[0]
        rpk = N / length
        return rpk / 1e6

    def rpkm(self, df_gene, df_rna):
        length = abs(df_gene['start'] - df_gene['end'])
        N = df_rna.shape[0]
        rpm = N / 1e6
        return rpm / length

    def check_empty(self):
        df = pd.read_csv('E-MTAB-1733.csv')
        fid = []
        for idx in df.index:
            url = df.loc[idx, 'link']
            _, fname = os.path.split(url)
            fname, _ = os.path.splitext(fname)
            fname = fname.split('_')[0]
            fid.append(fname)
        fid = np.unique(np.array(fid))

        dirname = os.path.join(self.root, 'database/RNA-seq')
        flist = os.listdir(dirname)
        flist = [os.path.splitext(x)[0] for x in flist if x.endswith('.db')]

        remains = list(set(fid) - set(flist))
        print(remains)

    def split_table(self, df):
        dfs = []
        df_chr = df.groupby('chromosome')
        for chr, df_sub in sorted(df_chr):
            for str, df_sub_sub in sorted(df_sub.groupby('strand')):
                dfs.append(df_sub_sub[['start', 'end']].values.flatten().astype(np.int32))
        return dfs

    def set_data(self, dfs):
        data_length = np.zeros(len(dfs)).astype(np.int32)
        for i, df in enumerate(dfs):
            data_length[i] = df.shape[0]
        dfs = np.concatenate(dfs)

        data_length_cum = np.zeros(data_length.shape[0] + 1)
        data_length_cum[1:] = data_length.cumsum()
        return dfs, data_length_cum

    def to_gpu(self, dfs_ref, data_length_cum_ref, dfs_rna, data_length_cum_src, N):
        THREADS_PER_BLOCK = 1 << 10

        ref_buffer_gpu = cuda.mem_alloc(dfs_ref.nbytes)
        cuda.memcpy_htod(ref_buffer_gpu, dfs_ref)

        res_buffer_gpu = cuda.mem_alloc(dfs_rna.nbytes)
        cuda.memcpy_htod(res_buffer_gpu, dfs_rna)

        ref_data_lengths_cum_gpu = cuda.mem_alloc(data_length_cum_ref.nbytes)
        cuda.memcpy_htod(ref_data_lengths_cum_gpu, data_length_cum_ref)

        res_data_lengths_cum_gpu = cuda.mem_alloc(data_length_cum_src.nbytes)
        cuda.memcpy_htod(res_data_lengths_cum_gpu, data_length_cum_src)

        M = np.int32(data_length_cum_src.shape[0])

        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)
        out_buffer = np.zeros((N, 2)).flatten().astype(np.int32)
        out_buffer_gpu = cuda.mem_alloc(out_buffer.nbytes)

        func = mod.get_function("cuda_scanner")
        func(ref_buffer_gpu, ref_data_lengths_cum_gpu, res_buffer_gpu, res_data_lengths_cum_gpu, out_buffer_gpu, N, M,
             block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out_buffer, out_buffer_gpu)
        out_buffer = out_buffer.reshape((-1, 2))
        return pd.DataFrame(data=out_buffer, columns=['TABLE_NUM', 'RPKM'])

    def run(self):
        fpath_out = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.out.db')
        con_out = sqlite3.connect(fpath_out)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.tissues.db')
        con = sqlite3.connect(fpath)
        tissues = Database.load_tableList(con)

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq.db')
        con_rna = sqlite3.connect(fpath)

        for tissue in tissues:
            print(tissue)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con)
            if df.shape[1] <= 5:
                continue
            dfs = self.split_table(df)
            dfs_ref, data_length_cum_ref = self.set_data(dfs)

            df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con_rna)
            dfs_rna = self.split_table(df_rna)
            dfs_rna, data_length_cum_src = self.set_data(dfs_rna)

            df_out = self.to_gpu(dfs_ref, data_length_cum_ref, dfs_rna, data_length_cum_src, df.shape[0])
            df = pd.concat([df, df_out], axis=1)
            df.to_sql(tissue, con_out, if_exists='replace', index=None)


if __name__ == '__main__':
    fr = Fantom_RNA()
    fr.run()
