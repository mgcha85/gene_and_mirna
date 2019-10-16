import pandas as pd

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy as np


mod = SourceModule("""
#include <stdio.h>

enum{START, END, WIDTH};

__device__ int binary_search(const int *X, const int val, const int N, const int col, const int idx)
{
    // X: array, val: exact value, N: length of X, col: columns number
    // if X has no val, return -1 

    if(N <= 1) return -1;
    if(val <= X[col]) return 0;
    if(val >= X[WIDTH * (N - 1) + col]) return N;

    int start = 0;
    int end = N - 1;
    
    while(start <= end) 
    {
        int mid = (start + end) / 2;
        int res_val = X[mid * WIDTH + col];
        
        if(res_val == val)      return mid;
        else if(res_val < val)  start = mid + 1;
        else                    end = mid - 1;
    }
    return -1;
}

__device__ int binary_search_pos(const int *X, const int val, const int N, const int col, const int idx)
{
    // X: array, val: nth value, N: length of X, col: columns number
    // if X has no val, return -1 

    if(N <= 1) return -1;    
    if(val < X[START]) return -1;
    if(val > X[WIDTH * (N - 1) + END]) return -1;

    int start = 0;
    int end = N - 1;
    
    while(start < end) 
    {
        int mid = (start + end) / 2;
        int res_val = X[mid * WIDTH + col];
        
        if(res_val == val)      return mid;
        else if(res_val < val)  start = mid + 1;
        else                    end = mid - 1;
        
        if(start >= end) return start;
    }
    return -1;
}

__device__ void is_overlap(const int *X, const int ref_start, const int ref_end, int *out, int idx)
{
    if(out[START] < 0 && out[END] < 0) return;
    if(out[START] < 0) out[START] = 0;
    
    int out_start = out[START];
    int out_end = out[END] - 1;
    
    if(idx < 10) printf("[%d] out_start: %d, out_end: %d, ref_start: %d, ref_end: %d\\n", idx, out_start, out_end, ref_start, ref_end);
    for(int i=out_start; i<=out_end; i++) {
        if(idx == 2) printf("[%d] out_start: %d, out_end: %d, ref_start: %d, ref_end: %d, res_start: %d, res_end: %d\\n", idx, out_start, out_end, ref_start, ref_end, X[WIDTH * i + START], X[WIDTH * i + END]);
        if(!(X[WIDTH * i + START] > ref_end) && !(X[WIDTH * i + END] < ref_start))
            continue;
        else {
            out[START] = -1;
            out[END] = -1;
        }  
    }
}

__global__ void cuda_sql(const int *ref, const int *res, int *out, const int N, const int M)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N) return;
    
    int ref_start = ref[idx * WIDTH + START];
    int ref_end = ref[idx * WIDTH + END];
    
    for(int i=START; i<=END; i++)
        out[WIDTH * idx + i] = binary_search_pos(res, ref[idx * WIDTH + i], M, START, idx);
    
    is_overlap(res, ref_start, ref_end, &out[WIDTH * idx], idx);
}

""")


class Sqlgpu:
    def cpu_run(self, df_ref, df_res):
        out = pd.DataFrame(data=-np.ones((df_ref.shape[0], 2)), index=df_ref.index, columns=['start', 'end'])
        for idx in df_ref.index:
            if idx > 10:
                break
            start = df_ref.loc[idx, 'start']
            end = df_ref.loc[idx, 'end']
            df_res_sub = df_res[~(df_res['start'] > end) & ~(df_res['end'] < start)]
            if df_res_sub.empty:
                continue
            out.loc[idx, 'start'] = df_res_sub.index[0]
            out.loc[idx, 'end'] = df_res_sub.index[-1]
        return out

    def run(self, df_ref, df_res):
        # df_ref, df_res should be DataFrame
        # they should have start, end only
        # they should be separated by chromosome and strand
        # they should be sorted already

        THREADS_PER_BLOCK = 1 << 10

        # df_ref: [start, end]
        ref = df_ref.values.flatten().astype(np.int32)
        ref_gpu = cuda.mem_alloc(ref.nbytes)
        cuda.memcpy_htod(ref_gpu, ref)

        # df_res: [start, end, score]
        res = df_res.values.flatten().astype(np.int32)
        res_gpu = cuda.mem_alloc(res.nbytes)
        cuda.memcpy_htod(res_gpu, res)

        N = np.int32(df_ref.shape[0])
        M = np.int32(df_res.shape[0])

        # output
        out = -np.ones((N, 2)).flatten().astype(np.int32)
        out_gpu = cuda.mem_alloc(out.nbytes)

        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)

        func = mod.get_function("cuda_sql")
        func(ref_gpu, res_gpu, out_gpu, N, M, block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out, out_gpu)
        return np.sort(out.reshape((N, 2)), axis=1)


if __name__ == '__main__':
    import socket, os, sqlite3
    from Database import Database

    hostname = socket.gethostname()
    if hostname == 'mingyu-Precision-Tower-7810':
        root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
    elif hostname == 'DESKTOP-DLOOJR6':
        root = 'D:/Bioinformatics'
    else:
        root = '/lustre/fs0/home/mcha/Bioinformatics'

    # reference GENCODE
    ref_path = os.path.join(root, 'database/gencode', 'gencode.v30lift37.annotation_spt.db')
    con_ref = sqlite3.connect(ref_path)

    # Fantom5
    fan_path = os.path.join(root, 'database/Fantom/v5/tissues', 'fantom_cage_by_tissue_spt.db')
    con_fan = sqlite3.connect(fan_path)

    # RNA-seq
    rna_path = os.path.join(root, 'database/RNA-seq/out', 'RNA_seq_tissue_spt.db')
    con_rna = sqlite3.connect(rna_path)

    tlist = Database.load_tableList(con_ref)
    sqlgpu = Sqlgpu()
    for tname in tlist[:1]:
        source, version, chromosome, strand = tname.split('.')
        df_ref = pd.read_sql_query("SELECT start, end FROM '{}'".format(tname), con_ref)
        df_res = pd.read_sql_query("SELECT start, end FROM 'adrenal_{}_{}'".format(chromosome, strand), con_rna)

        # out = sqlgpu.cpu_run(df_ref, df_res)
        out = sqlgpu.run(df_ref, df_res)
        pd.DataFrame(out, columns=['start', 'end']).to_excel(sqlgpu.__class__.__name__ + '.xlsx', index=None)
