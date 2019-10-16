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

__device__ int binary_search_loc(const int *X, const int val, const int N, const int col, const int idx)
{
    // X: array, val: nth value, N: length of X, col: columns number
    // if X has no val, return -1 

    if(N <= 1) return -1;    
    if(val < X[col]) return -1;
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
        
        if(start >= end) return start-1;
    }
    return -1;
}

__device__ void get_overlap(const int *X, const int ref_start, const int ref_end, int *addr, const int N, const int idx)
{
    int offset = -1;
    int num_ele = 0;
    
    if((X[START] > ref_end) || (X[WIDTH * (N - 1) + END] < ref_start)) return;
    
    for(int i=0; i<N; i++) {
        int res_start = X[WIDTH * i + START];
        int res_end = X[WIDTH * i + END];
        
        if(!(res_start > ref_end) && !(res_end < ref_start)) {
            if(offset < 0)
                offset = i;
            num_ele++;        
        }
        if(ref_end < res_start) break;
    }
    if(offset < 0) return;
    
    addr[WIDTH * idx + START] = offset;
    addr[WIDTH * idx + END] = offset + num_ele - 1;    
}

__device__ void cuda_sum(const float *score, const int *addr, float *out, const int N, const int idx)
{
    for(int i=0; i<N; i++) {
        int out_start = addr[WIDTH * i + START];   
        int out_end = addr[WIDTH * i + END];   
        if(out_start < 0) continue;
        
        float sum = 0.0;
        for(int j=out_start; j<=out_end; j++) {
            sum += score[j];
        }
        if(sum > 0) out[i] = sum;
    }
}

__global__ void cuda_sql(const int *ref, const int *res, int *addr, float *score, float *out, const int N, const int M)
{
    const int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N) return;
    
    int ref_start = ref[idx * WIDTH + START];
    int ref_end = ref[idx * WIDTH + END];
    
    get_overlap(res, ref_start, ref_end, addr, M, idx);
    cuda_sum(score, addr, out, N, idx);
}

""")


class Sqlgpu:
    def cpu_run(self, df_ref, df_res):
        out = pd.DataFrame(data=-np.ones((df_ref.shape[0], 2)), index=df_ref.index, columns=['start', 'end'])
        for idx in df_ref.index:
            if idx > 100:
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
        res = df_res[['start', 'end']].values.flatten().astype(np.int32)
        res_gpu = cuda.mem_alloc(res.nbytes)
        cuda.memcpy_htod(res_gpu, res)

        score = df_res[['score']].values.flatten().astype(np.float32)
        score_gpu = cuda.mem_alloc(score.nbytes)
        cuda.memcpy_htod(score_gpu, score)

        N = np.int32(df_ref.shape[0])
        M = np.int32(df_res.shape[0])

        # output
        addr = (-1 * np.ones(N * 2)).astype(np.int32)
        addr_gpu = cuda.mem_alloc(addr.nbytes)
        cuda.memcpy_htod(addr_gpu, addr)

        out = np.zeros(N).astype(np.float32)
        out_gpu = cuda.mem_alloc(out.nbytes)
        cuda.memcpy_htod(out_gpu, out)

        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)

        func = mod.get_function("cuda_sql")
        func(ref_gpu, res_gpu, addr_gpu, score_gpu, out_gpu, N, M, block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(addr, addr_gpu)
        cuda.memcpy_dtoh(out, out_gpu)
        return addr.reshape((N, 2)), out


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
        print(tname)
        source, version, chromosome, strand = tname.split('.')
        df_ref = pd.read_sql_query("SELECT start, end FROM '{}'".format(tname), con_ref)
        df_res = pd.read_sql_query("SELECT start, end, FPKM FROM 'adrenal_{}_{}'".format(chromosome, strand), con_rna)
        df_res = df_res.rename(columns={"FPKM": "score"})

        # out = sqlgpu.cpu_run(df_ref, df_res)
        addr, out = sqlgpu.run(df_ref, df_res)
        pd.DataFrame(addr, columns=['start', 'end']).to_excel(sqlgpu.__class__.__name__ + '_addr.xlsx', index=None)
        pd.DataFrame(out, columns=['score']).to_excel(sqlgpu.__class__.__name__ + '_score.xlsx', index=None)
