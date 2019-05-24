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
enum{TABLE_NUM, FAN_START, FAN_END, OUT_WIDTH};
enum{START, END, WIDTH};

__device__ int* search(int *ref_buffer, int *location, const int ref_tss, const int idx, const int N)
{
    int ref_start = ref_tss - 100;
    int ref_end = ref_tss + 100;
    int res_start, res_end;
    
    for(int i=0; i<N; i++) {
        res_start = ref_buffer[i * WIDTH + START];
        res_end = ref_buffer[i * WIDTH + END];
        if(ref_start > res_end || ref_end < res_start)  continue;
        else {
            location[0] = res_start;
            location[1] = res_end;
            return location;
        }
                
    }
    return location;
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
    
    int location[2] = {0, };
    search(&res_buffer_gpu[sidx * WIDTH], location, tss, idx, eidx - sidx);
    out_buffer_gpu[idx * OUT_WIDTH + FAN_START] = location[0];
    out_buffer_gpu[idx * OUT_WIDTH + FAN_END] = location[1];
    
}""")


class Comparison:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.chr_str_map = {}

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

    def set_ref_data(self, con):
        tlist = Database.load_tableList(con)

        buffer = []
        data_lengths = []
        for i, tname in enumerate(sorted(tlist)):
            _, chromosome, strand = tname.split('_')
            self.chr_str_map['{}_{}'.format(chromosome, strand)] = i
            df = pd.read_sql_query("SELECT start, end FROM '{}'".format(tname), con)
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

    def set_data(self, df_res):
        df_chr = df_res.groupby('chromosome')

        buffer = [None] * len(self.chr_str_map)
        data_lengths = [None] * len(self.chr_str_map)
        for chr, df_str in df_chr:
            for str, df_sub in df_str:
                num = self.chr_str_map['{}_{}'.format(chr, str)]
                buffer[num] = df_sub[['start', 'end']].values.flatten().astype(np.int32)
                data_lengths[num] = df_sub.shape[0]

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
        out_buffer = np.zeros((N, 3)).flatten().astype(np.int32)
        out_buffer_gpu = cuda.mem_alloc(out_buffer.nbytes)

        func = mod.get_function("cuda_scanner")
        func(ref_buffer_gpu, ref_data_lengths_cum_gpu, res_buffer_gpu, res_data_lengths_cum_gpu, out_buffer_gpu, N, M,
             block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out_buffer, out_buffer_gpu)
        out_buffer = out_buffer.reshape((N, 3))
        return out_buffer

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

    def run(self):
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation2.db')
        con_ref = sqlite3.connect(ref_path)
        ref_buffer, ref_data_lengths_cum = self.set_ref_data(con_ref)

        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue.db')
        con_fan = sqlite3.connect(fpath_fan)
        tlist_fan = Database.load_tableList(con_fan)

        for tname in tlist_fan:
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_fan)
            res_buffer, res_data_lengths_cum = self.set_data(df_fan)
            out_data = self.to_gpu(ref_buffer, ref_data_lengths_cum, res_buffer, res_data_lengths_cum)
            self.to_output(con_ref, out_data, ref_data_lengths_cum)


if __name__ == '__main__':
    comp = Comparison()
    if comp.hostname == 'mingyu-Precision-Tower-7810':
        # comp.to_server()
        comp.merge_table()
    else:
        comp.run()
