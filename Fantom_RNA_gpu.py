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

__device__ int get_table_num(int *ref_data_lengths_cum_gpu, const int idx, const int N)
{
    for(int i=0; i<N-1; i++) {  
        if(ref_data_lengths_cum_gpu[i] <= idx && ref_data_lengths_cum_gpu[i+1] > idx) return i;
    }
    return N;
}

__global__ void cuda_scanner(int *ref_buffer_gpu, int *ref_data_lengths_cum_gpu, int *res_buffer_gpu, int *res_data_lengths_cum_gpu, float *out_buffer_gpu, const int N, const int M)
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
    out_buffer_gpu[idx * OUT_WIDTH + RPKM] = get_rpkm(&res_buffer_gpu[sidx * WIDTH], ref_start, ref_end, eidx - sidx, idx);

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
        df_sorted = []
        dfs = []
        data_length = []
        df_chr = df.groupby('chromosome')
        chromosomes = sorted(df_chr.groups)
        for chr in chromosomes:
            if len(chr) > 5:
                continue
            df_sub = df_chr.get_group(chr)
            df_sub_str = df_sub.groupby('strand')
            strands = sorted(df_sub_str.groups)
            for str in strands:
                df_sub_sub = df_sub_str.get_group(str)
                dfs.append(df_sub_sub[['start', 'end']].values.flatten().astype(np.int32))
                data_length.append(df_sub_sub.shape[0])
                df_sorted.append(df_sub_sub)
        return dfs, data_length, pd.concat(df_sorted).reset_index(drop=True)

    def extract_tags_by_tissue(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.tissues.db')
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_out.db'))

        tlist = Database.load_tableList(con)
        dfs = []
        for tname in tlist:
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            if df.empty:
                continue
            df.loc[:, 'tissue'] = tname
            dfs.append(df)

        df_con = pd.concat(dfs)
        df_con['location'] = df_con['chromosome'] + ';' + df_con['start'].astype('str') + ';' + df_con['end'].astype('str') + ';' + df_con['strand']
        df_con_grp = df_con.groupby('location')

        result = []
        N = len(df_con_grp)
        for i, (group, df_sub) in enumerate(df_con_grp):
            if i % 1000 == 0 or i + 1 == N:
                print('{:0.2f}%'.format(100 * (i + 1) / N))

            tissues = ';'.join(df_sub['tissue'])
            result.append([*group.split(';'), tissues])

        df_res = pd.DataFrame(data=result, columns=['chromosome', 'start', 'end', 'strand', 'tissues'])
        df_res_chr = df_res.groupby('chromosome')
        for chr, df_sub in df_res_chr:
            df_sub_str = df_sub.groupby('strand')
            for str, df_sub_sub in df_sub_str:
                df_sub_sub = df_sub_sub.drop(['chromosome', 'strand'], axis=1)
                df_sub_sub[['start', 'end']] = df_sub_sub[['start', 'end']].astype(int)
                df_sub_sub.sort_values('start').to_sql('FANTOM_{}_{}'.format(chr, str), con_out, if_exists='replace', index=None)

    def split_by_tissue(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.fantom_cross_check_ucsc_ensembl_gencode.db')
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_out.db'))

        fpath_tis = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.tissues_out.db')
        con_tis = sqlite3.connect(fpath_tis)

        tlist = Database.load_tableList(con)
        for tname in tlist:
            print(tname)
            _, chromosome, strand = tname.split('_')
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            df_tis = pd.read_sql_query("SELECT * FROM 'FANTOM_{}_{}'".format(chromosome, strand), con_tis)
            df_tis.index = df_tis['start'].astype(str) + ';' + df_tis['end'].astype(str)

            tissues = []
            for idx in df.index:
                start = df.loc[idx, 'fan_start']
                end = df.loc[idx, 'fan_end']
                key = '{};{}'.format(start, end)
                if key in df_tis.index:
                    tissues.append(df_tis.loc[key, 'tissues'])
                else:
                    tissues.append('')

            df.loc[:, 'tissues'] = tissues
            df.to_sql(tname, con_out, if_exists='replace', index=None)

    def get_confirmed_fantom(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.fantom_cross_check_ucsc_ensembl_gencode.db')
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_out.db'))
        tlist = Database.load_tableList(con)

        for tname in tlist:
            print(tname)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            df_grp = df.groupby('gene_name')

            contents = []
            for gene, df_sub in df_grp:
                resources = np.unique(df_sub['resource'].values)
                df_sub = df_sub.drop_duplicates(subset=['fan_start', 'fan_end'])
                df_sub = df_sub[['fan_start', 'fan_end', 'tissues', 'gene_name']]
                df_sub.loc[:, 'resources'] = ';'.join(resources)
                contents.append(df_sub)
            df_res = pd.concat(contents)
            df_res.to_sql(tname, con_out, if_exists='replace', index=None)

    def split_by_tissues(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.fantom_cross_check_ucsc_ensembl_gencode_out.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        for tname in tlist:
            print(tname)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)

            result = {}
            for idx in df.index:
                tissues = df.loc[idx, 'tissues'].split(';')
                start = df.loc[idx, 'fan_start']
                end = df.loc[idx, 'fan_end']
                resources = df.loc[idx, 'resources']

                for tis in tissues:
                    if tis not in result:
                        result[tis] = [[start, end, resources]]
                    else:
                        result[tis].append([start, end, resources])

            for tis, contents in result.items():
                con_out = sqlite3.connect(fpath.replace('.db', '_{}.db'.format(tis)))
                df_res = pd.DataFrame(data=contents, columns=['start', 'end', 'resources'])
                df_res.to_sql(tname, con_out, if_exists='replace', index=None)

    def set_data(self, dfs, data_length):
        dfs = np.concatenate(dfs)
        data_length = np.array(data_length)
        data_length_cum = np.zeros(data_length.shape[0] + 1)
        data_length_cum[1:] = data_length.cumsum()
        return dfs, data_length_cum.astype(np.int32)

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

        M = np.int32(data_length_cum_ref.shape[0])

        gridN = int((N + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK)
        out_buffer = np.zeros((N, 2)).flatten().astype(np.float32)
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
        tissues_src = Database.load_tableList(con_rna)

        for tissue in tissues:
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con)
            if df.shape[1] <= 5:
                continue
            dfs, data_length, df_sorted = self.split_table(df)
            dfs_ref, data_length_cum_ref = self.set_data(dfs, data_length)

            tsrc = [x for x in tissues_src if tissue in x]
            for t in tsrc:
                df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(t), con_rna)
                df_rna = df_rna.rename(columns={'stop': 'end'})
                df_rna = df_rna[df_rna['strand'] != '.']
                dfs_rna, data_length_rna, df_sorted_rna = self.split_table(df_rna)
                dfs_rna, data_length_cum_src = self.set_data(dfs_rna, data_length_rna)

                df_out = self.to_gpu(dfs_ref, data_length_cum_ref, dfs_rna, data_length_cum_src, np.int32(df.shape[0]))
                df_out = pd.concat([df_sorted, df_out], axis=1)
                df_out.to_sql(t, con_out, if_exists='replace', index=None)
                print('{} is done'.format(t))


if __name__ == '__main__':
    fr = Fantom_RNA()
    fr.get_confirmed_fantom()
    fr.split_by_tissues()
