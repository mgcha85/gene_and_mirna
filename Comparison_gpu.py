import sqlite3
import pandas as pd
import socket
import os
import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
from Database import Database

mod = SourceModule("""
#include <stdio.h>
enum{TABLE_NUM, RESULT, OUT_WIDTH};
enum{START, END, WIDTH};

__device__ int search(int *ref_buffer, const int ref_tss, const int idx, const int N)
{
    int ref_start = ref_tss - 100;
    int ref_end = ref_tss + 100;
    int res_start, res_end;
    
    for(int i=0; i<N; i++) {
        res_start = ref_buffer[i * WIDTH + START];
        res_end = ref_buffer[i * WIDTH + END];
        if(ref_start > res_end || ref_end < res_start)  continue;
        else    return 1;        
    }
    return 0;
}

__device__ int get_table_num(int *ref_data_lengths_cum_gpu, const int idx, const int N)
{
    for(int i=0; i<N-1; i++) {  
        if(ref_data_lengths_cum_gpu[i] <= idx && ref_data_lengths_cum_gpu[i+1] > idx) return i;
    }
    return N;
}

__global__ void cuda_scanner(int *ref_buffer_gpu, int *ref_data_lengths_cum_gpu, int *res_buffer_gpu, int *res_data_lengths_cum_gpu, int *out_buffer_gpu, const int N, const int M)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N) return;

    int tb_num, sidx, eidx, tss;

    tss = ref_buffer_gpu[idx];
    tb_num = get_table_num(ref_data_lengths_cum_gpu, idx, M);
    sidx = res_data_lengths_cum_gpu[tb_num];
    eidx = res_data_lengths_cum_gpu[tb_num + 1];
    
    out_buffer_gpu[idx * OUT_WIDTH + TABLE_NUM] = tb_num;
    out_buffer_gpu[idx * OUT_WIDTH + RESULT] = search(&res_buffer_gpu[sidx * WIDTH], tss, idx, eidx - sidx);
    
}""")


class Comparison:
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
        self.chr_str_map = {}

    def set_ref_data(self, con):
        tlist = Database.load_tableList(con)

        buffer = []
        data_lengths = []
        for i, tname in enumerate(sorted(tlist)):
            _, chromosome, strand = tname.split('_')
            self.chr_str_map['{}_{}'.format(chromosome, strand)] = i
            df = pd.read_sql_query("SELECT start, end FROM '{}'".format(tname), con)
            # df = pd.read_sql_query("SELECT start, end FROM '{}' WHERE feature='transcript'".format(tname), con)
            if strand == '+':
                df['tss'] = df['start']
            else:
                df['tss'] = df['end']
            buffer.append(df['tss'].values.astype(np.int32))
            data_lengths.append(df.shape[0])

        data_lengths = np.array(data_lengths)
        data_lengths_cum = np.zeros(data_lengths.shape[0] + 1)
        data_lengths_cum[1:] = data_lengths.cumsum()
        return buffer, data_lengths_cum.astype(np.int32)

    def set_data(self, con):
        tlist = Database.load_tableList(con)

        buffer = []
        data_lengths = []
        for i, tname in enumerate(sorted(tlist)):
            _, chromosome, strand = tname.split('_')
            df = pd.read_sql_query("SELECT start, end FROM '{}'".format(tname), con)
            buffer.append(df.values.flatten().astype(np.int32))
            data_lengths.append(df.shape[0])

        data_lengths = np.array(data_lengths)
        data_lengths_cum = np.zeros(data_lengths.shape[0] + 1)
        data_lengths_cum[1:] = data_lengths.cumsum()
        return buffer, data_lengths_cum.astype(np.int32)

    def to_gpu(self, ref_buffer, ref_data_lengths_cum, res_buffer, res_data_lengths_cum):
        THREADS_PER_BLOCK = 1 << 10

        ref_buffer = np.concatenate(ref_buffer)
        res_buffer = np.concatenate(res_buffer)

        ref_buffer_gpu = cuda.mem_alloc(ref_buffer.nbytes)
        cuda.memcpy_htod(ref_buffer_gpu, ref_buffer)

        res_buffer_gpu = cuda.mem_alloc(res_buffer.nbytes)
        cuda.memcpy_htod(res_buffer_gpu, res_buffer)

        ref_data_lengths_cum_gpu = cuda.mem_alloc(ref_data_lengths_cum.nbytes)
        cuda.memcpy_htod(ref_data_lengths_cum_gpu, ref_data_lengths_cum)

        res_data_lengths_cum_gpu = cuda.mem_alloc(res_data_lengths_cum.nbytes)
        cuda.memcpy_htod(res_data_lengths_cum_gpu, res_data_lengths_cum)

        N = np.int32(ref_buffer.shape[0])
        M = np.int32(ref_data_lengths_cum.shape[0])

        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)
        out_buffer = np.zeros((N, 2)).flatten().astype(np.int32)
        out_buffer_gpu = cuda.mem_alloc(out_buffer.nbytes)

        func = mod.get_function("cuda_scanner")
        func(ref_buffer_gpu, ref_data_lengths_cum_gpu, res_buffer_gpu, res_data_lengths_cum_gpu, out_buffer_gpu, N, M,
             block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out_buffer, out_buffer_gpu)
        out_buffer = out_buffer.reshape((N, 2))
        return out_buffer

    def to_output(self, con_ref, out_data, data_lengths_cum, key):
        tlist = Database.load_tableList(con_ref)
        for i, tname in enumerate(sorted(tlist)):
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_ref)
            sidx = data_lengths_cum[i]
            eidx = data_lengths_cum[i + 1]

            df['table_num'] = out_data[sidx: eidx, 0]
            df[key] = out_data[sidx: eidx, 1]
            df.to_sql(tname, con_ref, if_exists='replace', index=None)

    def extract_attribute(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation.db')
        con = sqlite3.connect(fpath)

        fpath_out = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation_pc.db')
        con_out = sqlite3.connect(fpath_out)

        tlist = Database.load_tableList(con)
        index = []
        additional = []
        for tname in tlist:
            print(tname)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            attribute = df['attribute'].str.split('; ')
            for idx in attribute.index:
                attr = attribute.loc[idx]
                gene_id = attr[0].split(' ')[1].replace('"', '')
                transcript_id = attr[1].split(' ')[1].replace('"', '')
                gene_name = attr[3].split(' ')[1].replace('"', '')
                transcript_type = attr[4].split(' ')[1].replace('"', '')
                transcript_name = attr[5].split(' ')[1].replace('"', '')
                if transcript_type == 'protein_coding':
                    index.append(idx)
                    additional.append([gene_id, transcript_id, gene_name, transcript_name])

            df_add = pd.DataFrame(data=additional, columns=['gene_id', 'transcript_id', 'gene_name', 'transcript_name'])
            df = df.loc[index].reset_index(drop=True)
            df = pd.concat([df, df_add], axis=1)
            df.to_sql(tname, con_out, if_exists='replace', index=None)

    def run(self):
        ref_path = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.db')
        # ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation.db')
        con_ref = sqlite3.connect(ref_path)

        fpath_ens = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19.db')
        fpath_ucsc = os.path.join(self.root, 'database/UCSC/Genes', 'genes.db')
        fpath_fan = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation.db')

        con_ens = sqlite3.connect(fpath_ens)
        con_ucsc = sqlite3.connect(fpath_ucsc)
        con_fan = sqlite3.connect(fpath_fan)

        ref_buffer, ref_data_lengths_cum = self.set_ref_data(con_ref)
        for key, con_res in zip(['ENS', 'UCSC', 'GENCODE'], [con_ens, con_ucsc, con_fan]):
            res_buffer, res_data_lengths_cum = self.set_data(con_res)
            print(key, res_data_lengths_cum[-1])
            out_data = self.to_gpu(ref_buffer, ref_data_lengths_cum, res_buffer, res_data_lengths_cum)
            self.to_output(con_ref, out_data, ref_data_lengths_cum, key)

    def stats(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.db')
        # fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        dfs = []
        for tname in tlist:
            dfs.append(pd.read_sql_query("SELECT * FROM '{}'".format(tname), con))
        df_res = pd.concat(dfs)

        stats = {'ENS': 0, 'UCSC': 0, 'GENCODE': 0}
        tot = df_res.shape[0]
        for col, _ in stats.items():
            npos = df_res[df_res[col] > 0].shape[0]
            print(npos, tot)
            stats[col] = 100 * npos / tot
        print(stats)


if __name__ == '__main__':
    comp = Comparison()
    comp.extract_attribute()
    # comp.stats()
