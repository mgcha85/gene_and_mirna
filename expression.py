import sqlite3
import pandas as pd
import socket
import os
import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
from Database import Database
from Server import Server
import sys


mod = SourceModule("""
#include <stdio.h>
enum{LABEL, INDEX, OUT_WIDTH};
enum{START, END, WIDTH};
#define MAXLEN  1 << 10


__device__ int* search(int *res_buffer, int *out_buffer_gpu, const int ref_tss, const int idx, const int N, const int scope)
{
    int ref_start = ref_tss - scope;
    int ref_end = ref_tss + scope;
    int res_start, res_end;

    for(int i=0; i<N; i++) {
        res_start = res_buffer[i * WIDTH + START];
        res_end = res_buffer[i * WIDTH + END];
        if(ref_start > res_end || ref_end < res_start)  continue;
        else {
            out_buffer_gpu[OUT_WIDTH * i + LABEL] = 1;
            out_buffer_gpu[OUT_WIDTH * i + INDEX] = idx;
        }
    }
}

__device__ int get_table_num(int *ref_data_lengths_cum_gpu, const int idx, const int N)
{
    for(int i=0; i<N-1; i++) {
        if(ref_data_lengths_cum_gpu[i] <= idx && ref_data_lengths_cum_gpu[i+1] > idx) return i;
    }
    return N;
}

__global__ void cuda_scanner(int *ref_buffer_gpu, int *ref_data_lengths_cum_gpu, int *res_buffer_gpu, int *res_data_lengths_cum_gpu, int *out_buffer_gpu, const int N, const int M, const int scope)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N) return;

    int tb_num, sidx, eidx, tss;

    tss = ref_buffer_gpu[idx];
    tb_num = get_table_num(ref_data_lengths_cum_gpu, idx, M);
    sidx = res_data_lengths_cum_gpu[tb_num];
    eidx = res_data_lengths_cum_gpu[tb_num + 1];

    search(&res_buffer_gpu[sidx * WIDTH], &out_buffer_gpu[sidx * OUT_WIDTH], tss, idx, eidx - sidx, scope);

}""")


class Expression:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinfomatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.chr_str_map = {}

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

    def set_ref_data(self, con):
        tlist = Database.load_tableList(con)

        dfs = []
        buffer = []
        data_lengths = []
        for i, tname in enumerate(sorted(tlist)):
            chromosome, strand = tname.split('_')[-2:]
            if chromosome == 'chrY':
                continue
            self.chr_str_map['{}_{}'.format(chromosome, strand)] = i
            df = pd.read_sql_query("SELECT start, end, gene_name FROM '{}'".format(tname), con)
            if df.empty:
                continue

            if strand == '+':
                df['tss'] = df['start']
            else:
                df['tss'] = df['end']
            buffer.append(df['tss'].values.astype(np.int32))
            data_lengths.append(df.shape[0])
            dfs.append(df)

        data_lengths = np.array(data_lengths)
        data_lengths_cum = np.zeros(data_lengths.shape[0] + 1)
        data_lengths_cum[1:] = data_lengths.cumsum()
        return buffer, data_lengths_cum.astype(np.int32), pd.concat(dfs).reset_index(drop=True)

    def set_data(self, df_res):
        df_chr = df_res.groupby('chromosome')

        N = len(self.chr_str_map)
        dfs = [None] * N
        buffer = [None] * N
        data_lengths = [None] * N
        for chr_str in self.chr_str_map.keys():
            chr, str = chr_str.split('_')
            df_str = df_chr.get_group(chr)
            for str, df_sub in df_str.groupby('strand'):
                key = '{}_{}'.format(chr, str)
                if key not in self.chr_str_map:
                    continue

                if chr == 'chrMT':
                    chr = 'chrM'
                if chr == 'chrY' or str == '.' or len(chr) > 5:
                    continue
                num = self.chr_str_map['{}_{}'.format(chr, str)]
                buffer[num] = df_sub[['start', 'end']].values.flatten().astype(np.int32)
                data_lengths[num] = df_sub.shape[0]
                dfs[num] = df_sub

        data_lengths = np.array(data_lengths)
        data_lengths_cum = np.zeros(data_lengths.shape[0] + 1)
        data_lengths_cum[1:] = data_lengths.cumsum()
        return buffer, data_lengths_cum.astype(np.int32), pd.concat(dfs).reset_index(drop=True)

    def to_gpu(self, ref_buffer, ref_data_lengths_cum, res_buffer, res_data_lengths_cum, scope):
        THREADS_PER_BLOCK = 1 << 10
        SCOPE = np.int32(scope)

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
        out_buffer = np.zeros((len(res_buffer) // 2, 2)).flatten().astype(np.int32)
        out_buffer_gpu = cuda.mem_alloc(out_buffer.nbytes)

        func = mod.get_function("cuda_scanner")
        func(ref_buffer_gpu, ref_data_lengths_cum_gpu, res_buffer_gpu, res_data_lengths_cum_gpu, out_buffer_gpu, N, M,
             SCOPE, block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out_buffer, out_buffer_gpu)

        df_res = pd.DataFrame(res_buffer.reshape((-1, 2)), columns=['start', 'end'])
        df_out = pd.DataFrame(out_buffer.reshape((-1, 2)), columns=['label', 'index'])
        return pd.concat([df_res, df_out], axis=1)

    def to_output(self, con_ref, out_data, data_lengths_cum):
        tlist = Database.load_tableList(con_ref)
        for i, tname in enumerate(sorted(tlist)):
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_ref)
            sidx = data_lengths_cum[i]
            eidx = data_lengths_cum[i + 1]

            df['table_num'] = out_data[sidx: eidx, 0]
            df['fan_start'] = out_data[sidx: eidx, 1]
            df['fan_end'] = out_data[sidx: eidx, 2]
            df.to_sql(tname, con_ref, if_exists='replace', index=None)

    def fantom_to_gene(self, scope=500):
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation2.db')
        con_ref = sqlite3.connect(ref_path)
        ref_buffer, ref_data_lengths_cum, df_ref = self.set_ref_data(con_ref)

        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE.db')
        con_fan = sqlite3.connect(fpath_fan)
        tlist_fan = Database.load_tableList(con_fan)

        fpath_fan_out = fpath_fan.replace('.db', '_{}.db'.format(scope))

        con_out = sqlite3.connect(fpath_fan_out)
        for i, tname in enumerate(tlist_fan):
            print(tname)
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_fan)
            res_buffer, res_data_lengths_cum, df_fan_sorted = self.set_data(df_fan)
            df_out = self.to_gpu(ref_buffer, ref_data_lengths_cum, res_buffer, res_data_lengths_cum, scope=scope)
            df_res = pd.concat([df_fan_sorted, df_out[['label', 'index']]], axis=1)
            df_res = df_res[df_res['label'] > 0].reset_index(drop=True)
            if df_res.empty:
                continue
            gene_name = df_ref.loc[df_res['index'], 'gene_name'].values
            df_res.loc[:, 'gene_name'] = gene_name
            df_res.drop(['label', 'index'], axis=1).to_sql(tname, con_out, if_exists='replace', index=None)

    def fantom_to_mir(self, scope=500):
        ref_path = os.path.join(self.root, 'database', 'fantom5_mir.db')
        con_ref = sqlite3.connect(ref_path)
        ref_buffer, ref_data_lengths_cum, df_ref = self.set_ref_data(con_ref)

        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE.db')
        con_fan = sqlite3.connect(fpath_fan)
        tlist_fan = Database.load_tableList(con_fan)
        con_out = sqlite3.connect(fpath_fan.replace('.db', '_mir_{}.db'.format(scope)))

        for i, tname in enumerate(tlist_fan):
            print(tname)
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_fan)
            res_buffer, res_data_lengths_cum, df_fan_sorted = self.set_data(df_fan)
            df_out = self.to_gpu(ref_buffer, ref_data_lengths_cum, res_buffer, res_data_lengths_cum, scope=scope)
            df_res = pd.concat([df_fan_sorted, df_out[['label', 'index']]], axis=1)
            df_res = df_res[df_res['label'] > 0].reset_index(drop=True) # if any tags exist

            gene_name = [None] * df_res.shape[0]
            for i, idx in zip(df_res.index, df_res['index']):
                gene_name[i] = df_ref.loc[idx, 'gene_name']
            df_res.loc[:, 'gene_name'] = gene_name
            df_res.drop(['label', 'index'], axis=1).to_sql(tname, con_out, if_exists='replace', index=None)

    def rna_to_gene(self):
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation2.db')
        con_ref = sqlite3.connect(ref_path)
        ref_buffer, ref_data_lengths_cum, df_ref = self.set_ref_data(con_ref)

        fpath_fan = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con_fan = sqlite3.connect(fpath_fan)
        tlist_fan = Database.load_tableList(con_fan)
        con_out = sqlite3.connect(fpath_fan.replace('.db', '_out.db'))

        for i, tname in enumerate(tlist_fan):
            print(tname)
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_fan)
            res_buffer, res_data_lengths_cum, df_fan_sorted = self.set_data(df_fan)
            df_out = self.to_gpu(ref_buffer, ref_data_lengths_cum, res_buffer, res_data_lengths_cum)
            df_res = pd.concat([df_fan_sorted, df_out[['label', 'index']]], axis=1)
            df_res = df_res[df_res['label'] > 0].reset_index(drop=True)

            gene_name = [None] * df_res.shape[0]
            for i, idx in zip(df_res.index, df_res['index']):
                gene_name[i] = df_ref.loc[idx, 'gene_name']
            df_res.loc[:, 'gene_name'] = gene_name
            df_res.drop(['label', 'index'], axis=1).to_sql(tname, con_out, if_exists='replace', index=None)
            # self.to_output(con_ref, out_data, ref_data_lengths_cum)


if __name__ == '__main__':
    comp = Expression()
    if comp.hostname == 'mingyu-Precision-Tower-7810':
        comp.to_server()
    else:

        for scope in [100, 500]:
            for func in [comp.fantom_to_gene, comp.fantom_to_mir]:
                func(scope)
