import socket
import sqlite3
import os
import pandas as pd
import numpy as np
from XmlHandler import XmlHandler
import pickle as pkl
import matplotlib.pyplot as plt

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

mod = SourceModule("""
#include <stdio.h>

enum{CHROM, CHROM_START, CHROM_END, CHROM_WIDTH};
enum{CHROMOSOM, REF_START, REF_END, TSTART, TEND, TABLE_WIDTH};
enum{ADDR_OFFSET, ADDR_SIZE, ADDR_WIDTH};

__device__ int filterInRange(const int bin_start, const int bin_end, const int start, const int end) {
    if(!(bin_start > end) && !(bin_end < start))
        return 0;
    else
        return 1;
}
    
__device__ int get_hist_width(int ref_start, int ref_end, const int BINSIZE)
{
    int length = ref_end - ref_start;
    return (length + BINSIZE - 1) / BINSIZE;
}

__device__ void histogram(int ref_start, int ref_end, int *data_scr, const int *addr_table, int *out_histogram, int *out_buffer_size_cum, const int idx, const int BINSIZE)
{
    int offset = addr_table[idx * ADDR_WIDTH + ADDR_OFFSET];
    int size = addr_table[idx * ADDR_WIDTH + ADDR_SIZE];
    int hist_width = out_buffer_size_cum[idx + 1] - out_buffer_size_cum[idx];
    
    for(int j=0; j<size; j++) {
        int src_start = data_scr[(offset+j) * CHROM_WIDTH + CHROM_START];
        int src_end = data_scr[(offset+j) * CHROM_WIDTH + CHROM_END];
        for(int i=0; i<hist_width; i++) {
            int bin_start = ref_start + (i * BINSIZE);
            int bin_end = bin_start + BINSIZE;
            if(filterInRange(bin_start, bin_end, src_start, src_end) == 0)
                out_histogram[out_buffer_size_cum[idx] + i] += 1;
        }
    }    
}

__device__ void first_filter(int ref_start, int ref_end, int *data_scr, int table_offset, int table_size, int *addr_table, int idx)
{
    int table_end = table_offset + table_size;
    int cnt = 0;

    for (int i=table_offset; i<table_end; i++) {
        int scr_start = data_scr[i * CHROM_WIDTH + CHROM_START];
        int scr_end = data_scr[i * CHROM_WIDTH + CHROM_END];
        if(filterInRange(ref_start, ref_end, scr_start, scr_end) == 0) {
            if(addr_table[idx * ADDR_WIDTH + ADDR_OFFSET] < 0)
                addr_table[idx * ADDR_WIDTH + ADDR_OFFSET] = i;
            addr_table[idx * ADDR_WIDTH + ADDR_SIZE] = ++cnt;
        }
    }
}

__global__ void cuda_histogram(int *data_ref, int *data_scr, int *charom_table_size, int *addr_table, int *out_histogram, int *out_buffer_size_cum, const int BINSIZE, const int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N) return;

    int chrom = data_ref[idx*CHROM_WIDTH+CHROM];
    int start = data_ref[idx*CHROM_WIDTH+CHROM_START];
    int end = data_ref[idx*CHROM_WIDTH+CHROM_END];
    int table_offset = charom_table_size[chrom];
    int table_size = charom_table_size[chrom+1] - table_offset;

    first_filter(start, end, data_scr, table_offset, table_size, addr_table, idx);
    histogram(start, end, data_scr, addr_table, out_histogram, out_buffer_size_cum, idx, BINSIZE);

}""")


class histogram_gpu:
    def __init__(self, fpath='histogram_param.xml'):
        hostname = socket.gethostname()
        self.resource = 'gpu'
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/mnt/28218d84-f0fc-48a6-ae9d-f08cdc4d3a25/bioinformatics'
            self.resource = 'cpu'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
            self.resource = 'gpu'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
            self.resource = 'gpu'
        elif 'evc' in hostname:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
            self.resource = 'gpu'
        else:
            print('wrong option')
            return
        self.user_param = XmlHandler.load_param(fpath)

    def multiprocessing_sql(self, df_ref, dbpath_src, src_tname, idx):
        chromosome = df_ref.loc[idx, 'chromosome']
        start = df_ref.loc[idx, 'start']
        end = df_ref.loc[idx, 'end']

        con = sqlite3.connect(dbpath_src)
        tname = src_tname.format(chromosome)
        return pd.read_sql_query("SELECT * FROM '{}' WHERE NOT start>{end} AND NOT end<{start}"
                                 "".format(tname, start=start, end=end), con)

    def extract_from_sql(self, df_ref, dbpath_src, src_tname):
        contents = []
        for idx in df_ref.index:
            if idx % 100 == 0 or idx + 1 == df_ref.shape[0]:
                print('{:0.2f}%'.format(100 * (idx + 1) / df_ref.shape[0]))
            contents.append(self.multiprocessing_sql(df_ref, dbpath_src, src_tname, idx))
        return contents

    def load_tableList(self, con):
        cursor = con.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        df_slist = []
        for table_name in tables:
            df_slist.append(table_name[0])
        return np.array(sorted(df_slist))
    
    def get_out_buffer_size(self, df_ref, binsize):
        out_buffer_size = np.zeros(df_ref.shape[0]).astype(np.int32)
        for idx in df_ref.index:
            start = df_ref.loc[idx, 'start']
            end = df_ref.loc[idx, 'end']
            out_buffer_size[idx] = (end - start + binsize) // binsize

        out_buffer_size_cum = np.zeros(out_buffer_size.shape[0] + 1).astype(np.int32)
        out_buffer_size_cum[1:] = out_buffer_size.cumsum()
        return out_buffer_size_cum

    def resize_output(self, out_histogram, out_buffer_size_cum):
        out = []
        N = len(out_buffer_size_cum)
        for i in range(1, N):
            sidx = out_buffer_size_cum[i - 1]
            eidx = out_buffer_size_cum[i]
            out.append(out_histogram[sidx: eidx])
        return out

    def histogram_gpu(self, df_ref, dbpath_src, src_tname):
        THREADS_PER_BLOCK = 1 << 9
        N = np.int32(df_ref.shape[0])
        binsize = np.int32(self.user_param['numeric']['bin'])
        out_buffer_size_cum = self.get_out_buffer_size(df_ref, binsize)

        con = sqlite3.connect(dbpath_src)
        table_list = [x for x in self.load_tableList(con) if src_tname in x]

        data_scr = []
        charom_table_size = [0]
        chrom_table = {}
        chrom_table_by_str = {}
        for i, tname in enumerate(sorted(table_list)):
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            chrom_table[i] = tname.split('_')[1]
            chrom_table_by_str[tname.split('_')[1]] = i
            df['chromosome'] = i
            data_scr.append(df[['chromosome', 'start', 'end']].values.flatten())
            charom_table_size.append(df.shape[0])

        df_ref.loc[:, 'chromosome_'] = [chrom_table_by_str[x] for x in df_ref['chromosome']]

        # inputs to gpu
        data_ref = df_ref[['chromosome_', 'start', 'end']].values.flatten().astype(np.int32)
        data_scr = np.concatenate(data_scr).astype(np.int32)
        charom_table_size = np.array(charom_table_size).cumsum().astype(np.int32)
        addr_table = -np.ones((df_ref.shape[0], 2)).flatten().astype(np.int32)

        # output from gpu
        out_histogram = np.zeros(out_buffer_size_cum[-1]).flatten().astype(np.int32)

        # malloc to gpu
        data_ref_gpu = cuda.mem_alloc(data_ref.nbytes)
        data_scr_gpu = cuda.mem_alloc(data_scr.nbytes)
        charom_table_size_gpu = cuda.mem_alloc(charom_table_size.nbytes)
        addr_table_gpu = cuda.mem_alloc(addr_table.nbytes)
        out_histogram_gpu = cuda.mem_alloc(out_histogram.nbytes)

        # copy to gpu
        cuda.memcpy_htod(data_ref_gpu, data_ref)
        cuda.memcpy_htod(data_scr_gpu, data_scr)
        cuda.memcpy_htod(charom_table_size_gpu, charom_table_size)
        cuda.memcpy_htod(addr_table_gpu, addr_table)
        cuda.memcpy_htod(out_histogram_gpu, out_histogram)

        # set output buffer to gpu
        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)

        func = mod.get_function("cuda_histogram")
        func(data_ref_gpu, data_scr_gpu, charom_table_size_gpu, addr_table_gpu, out_histogram_gpu, out_buffer_size_cum, binsize, N,
             block=(THREADS_PER_BLOCK, 1, 1), grid=(gridN, 1))
        cuda.memcpy_dtoh(out_histogram, out_histogram_gpu)
        out = self.resize_output(out_histogram, out_buffer_size_cum)
        with open('out_hist.cha', 'wb') as f:
            pkl.dump(out, f)

        # out_histogram = out_histogram.reshape((-1, histogram_width))
        # df_add = pd.DataFrame(data=out_histogram, columns=range(histogram_width))
        # df_res = pd.concat([df_ref[['chromosome', 'start', 'end']], df_add], axis=1)
        # return df_res

    def run(self, df_ref, dbpath_src, src_tname):
        # hbw = int(self.user_param['numeric']['bandwidth'] // 2)
        #
        # df_ref.loc[:, 'tss'] = df_ref[['start', 'end']].mean(axis=1).astype(np.int32)
        # df_ref.loc[:, 'start'] = df_ref['tss'] - hbw
        # df_ref.loc[:, 'end'] = df_ref['tss'] + hbw
        return self.histogram_gpu(df_ref, dbpath_src, src_tname)

    def draw_pred_result(self):
        root = '/home/mingyu/Bioinformatics'
        user_param = XmlHandler.load_param('user_param.xml')

        for cline in ['A549']:
            # in_path = os.path.join(root, 'Papers/Trans', 'posResultsCNN.txt')
            # pd.read_csv(in_path, colu)
            dbpath_ref = os.path.join(root, 'Papers/mircoTSS/cover/Training', 'mitss_pred_score.db')
            dbpath_src = os.path.join(root, 'database/bioinfo_{}.db'.format(cline))
            tnames = user_param['table-names'][cline]
            con = sqlite3.connect(dbpath_ref)

            for tname in ['mitss_pred0', 'mitss_pred_fbd0', 'mitss_pred_histone0', 'mitss_pred_sequence0',
                          'mitss_pred_with_seq0', 'mitss_pred_with_seq_merge0', 'putative_tss_region0']:
                print(cline, tname)

                try:
                    df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
                except Exception as e:
                    print(e)
                    continue
                N = len(tnames)

                for i, (peak, tname2) in enumerate(tnames.items()):
                    # remove not common chromosome
                    df_ref = df_ref[df_ref['chromosome'] != 'chrM']
                    df_ref = df_ref[df_ref['chromosome'] != 'chrY']
                    df_ref = df_ref.reset_index()

                    df_res = self.run(df_ref, dbpath_src, tname2)
                    xaxis = (np.arange(df_res.shape[1] - 3)) * 20 - 500

                    plt.subplot(N, 1, i + 1)
                    plt.bar(xaxis, df_res.iloc[:, 3:].mean(axis=0), width=18)
                    plt.grid()
                    if i == 0:
                        plt.title('N: {:,d}'.format(df_ref.shape[0]))
                    else:
                        plt.xlabel('Distance from TSS (bp)')
                        plt.ylabel('Normalized Frequency')
                plt.savefig(os.path.join(root, tname + '_' + cline + '.png'))
                plt.close()


if __name__ == '__main__':
    hgpu = histogram_gpu()
    hostname = socket.gethostname()

    if hostname == 'mingyu-Precision-Tower-7810':
        root = '/home/mingyu/Bioinformatics'
    elif hostname == 'DESKTOP-DLOOJR6':
        root = 'D:/Bioinformatics'
    elif hostname == 'mingyu-Inspiron-7559':
        root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
    elif 'evc' in hostname:
        root = '/lustre/fs0/home/mcha/Bioinformatics'
    else:
        print('wrong option')
        exit(1)

    user_param = XmlHandler.load_param('user_param.xml')
    peak_opt = user_param['peak-option']

    for cline in ['GM12878', 'HelaS3', 'HepG2', 'K562']:
        in_path = os.path.join(root, 'papers_tss.db')

        dbpath_src = os.path.join(root, 'database/bioinfo_{}.db'.format(cline))
        tnames = user_param['table-names'][cline]
        con = sqlite3.connect(in_path)

        plt.figure(figsize=(12, 8))
        for tname, lstyle, label in zip(['cell_specific_{}', 'mitss_pred_with_seq'], ['-', '--', ':'], ['Hua', 'D-mirT']):
            tname = tname.format(cline)
            try:
                df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            except Exception as e:
                print(e)
                continue

            df_ref['tss'] = df_ref[['start', 'end']].mean(axis=1)
            df_ref['start'] = df_ref['tss']
            df_ref['end'] = df_ref['tss']

            N = len(tnames)
            for i, (peak, tname2) in enumerate(tnames.items()):
                # remove not common chromosome
                df_ref = df_ref[(df_ref['chromosome'] != 'chrM') & (df_ref['chromosome'] != 'chrY')]
                df_ref = df_ref.reset_index(drop=True)

                df_res = hgpu.run(df_ref, dbpath_src, tname2)
                xaxis = (np.arange(df_res.shape[1] - 3)) * 20 - 500

                plt.subplot(N, 1, i+1)
                plt.plot(xaxis, df_res.iloc[:, 3:].mean(axis=0), linestyle=lstyle, label=label, linewidth=2)
                plt.ylim([0, 0.5])
                if tname == 'mitss_pred_with_seq':
                    plt.grid()
                plt.xlabel('Distance from TSS (bp)')
                plt.ylabel('Normalized Frequency')

        plt.legend()
        plt.savefig(os.path.join(root, tname + '_' + cline + '.png'))
        plt.close()
