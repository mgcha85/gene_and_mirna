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
from joblib import Parallel, delayed
import multiprocessing


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

    __device__ int binary_search(float *X, int val, int N)
    {
        int start = 0;
        int end = N - 1;

        while(start < end) {
            if(start == end) return start;

            int mid = (start + end) / 2;
            if(X[mid] > val)    start = mid + 1;
            else    end = mid - 1;
        }
        return -1;
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
        out[idx] = correlationCoefficient(&X1[offset], &X2[offset], NUM_TISSUE);        
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
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.table_names = {}

    def to_server(self):
        from Server import Server
        import sys

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

    def processInput(self, tissue, loc):
        # Fantom5
        fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'fantom_cage_by_tissue_spt.db')
        con_fan = sqlite3.connect(fan_path)

        # RNA-seq
        rna_path = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_spt.db')
        con_rna = sqlite3.connect(rna_path)

        chromosome, strand, start, end = loc
        tname_rsc = '_'.join([tissue, chromosome, strand])
        df_fan = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND strand='{}' AND NOT "
                                   "start>{end} AND NOT end<{start}".format(tname_rsc, chromosome, strand,
                                                                            start=start, end=end), con_fan)
        df_rna = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND strand='{}' AND NOT "
                                   "start>{end} AND NOT end<{start}".format(tname_rsc, chromosome, strand,
                                                                            start=start, end=end), con_rna)
        if df_fan.empty and df_rna.empty:
            return [0, 0]

        return [df_fan['score'].sum(), df_rna['FPKM'].sum()]

    def correlation_fan_rna_cpu(self):
        num_cores = multiprocessing.cpu_count()
        if num_cores > 6:
            num_cores = 6
        print('num_cores: {}'.format(num_cores))

        # tissue list
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'tissues_fantom_rna.xlsx')
        df_tis = pd.read_excel(fpath, sheet_name='Sheet1')
        tissues = df_tis['RNA-seq']

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
        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_cpu.db')
        con_out = sqlite3.connect(out_path)

        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_cpu.db')
        con_sum_out = sqlite3.connect(out_path)

        columns = ['chromosome', 'start', 'end', 'strand', 'gene_name', 'transcript_name', 'transcript_type', 'corr']
        tlist = Database.load_tableList(ref_con)
        for tname in tlist:
            source, version, chromosome, strand = tname.split('.')
            if chromosome == 'chrM' or chromosome == 'chrY':
                continue

            print(chromosome, strand)

            df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)
            if strand == '+':
                offset = 'start'
            else:
                offset = 'end'

            tss = deepcopy(df_ref[offset])
            df_ref['end'] = tss + 500
            df_ref['start'] = tss - 500

            init_data = np.zeros((df_ref.shape[0], len(tissues) + 2))
            df_fan_sum = pd.DataFrame(data=init_data, index=df_ref.index, columns=['start', 'end', *tissues])
            df_fan_sum[['start', 'end']] = df_ref[['start', 'end']]
            df_rna_sum = pd.DataFrame(data=init_data, index=df_ref.index, columns=['start', 'end', *tissues])
            df_rna_sum[['start', 'end']] = df_ref[['start', 'end']]

            df_res = pd.DataFrame(index=df_ref.index, columns=columns)
            for idx in df_ref.index[20:]:
                if idx % 100 == 0 or idx + 1 == df_ref.shape[0]:
                    print('{:,d} / {:,d}'.format(idx + 1, df_ref.shape[0]))

                start = df_ref.loc[idx, 'start']
                end = df_ref.loc[idx, 'end']
                gene_name = df_ref.loc[idx, 'gene_name']
                transcript_name = df_ref.loc[idx, 'transcript_name']
                transcript_type = df_ref.loc[idx, 'transcript_type']

                df_buf = pd.DataFrame(data=np.zeros((M, 2)), index=tissues, columns=['FANTOM', 'RNA-seq'])
                for i, tissue in enumerate(tissues):
                    tname_rsc = '_'.join([tissue, chromosome, strand])
                    sql = "SELECT SUM(score) FROM '{}' WHERE chromosome='{} AND strand='{}' AND start<={end} AND NOT " \
                          "end>={start}".format(tname_rsc, chromosome, strand, start=start, end=end)
                    df_buf.loc[tissue, 'FANTOM'] = pd.read_sql_query(sql, con_fan)
                    df_buf.loc[tissue, 'RNA-seq'] = pd.read_sql_query(sql, con_rna)
                corr = np.corrcoef(df_buf['FANTOM'], df_buf['RNA-seq'])
                df_res.loc[idx, :] = [chromosome, start, end, strand, gene_name, transcript_name, transcript_type,
                                      corr[0, 1]]

                df_fan_sum.to_sql('_'.join(['fantom', chromosome, strand]), con_sum_out, index=None, if_exists='replace')
                df_rna_sum.to_sql('_'.join(['rna', chromosome, strand]), con_sum_out, index=None, if_exists='replace')
            df_res.to_sql(tname, con_out, if_exists='replace', index=None)

    def get_score_sum(self, df_ref, df_res):
        from sqlgpu import Sqlgpu

        sqlgpu = Sqlgpu()
        out = sqlgpu.bin_run(df_ref, df_res)
        return out

    def processInput_corr(self, df_ref, loc, cline, pct):
        print('{:.2f} %'.format(pct))
        fan_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_hCAGE_celllines.db')
        con_fan = sqlite3.connect(fan_path)
        chromosome, strand = loc

        contents = np.zeros(df_ref.shape[0])
        for i, idx in enumerate(df_ref.index):
            start = df_ref.loc[idx, 'start']
            end = df_ref.loc[idx, 'end']
            sql = "SELECT chromosome, strand, start, end, score FROM '{}' WHERE chromosome='{}' AND strand='{}' AND " \
                  "start<={} AND end>={}".format(cline, chromosome, strand, end, start)
            df = pd.read_sql_query(sql, con_fan).sort_values(by=['start'])
            contents[i] = df['score'].sum()
        return contents

    def correlation_fan_cpu(self, hbw, ref='gene'):
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count() - 2

        if ref == 'mir':
            # reference miRNA
            ref_path = os.path.join(self.root, 'database', 'consistent_miRNA_330.db')
            ref_con = sqlite3.connect(ref_path)
        else:
            # reference GENCODE
            ref_path = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
            ref_con = sqlite3.connect(ref_path)

        # Fantom5
        fan_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_hCAGE_celllines.db')
        con_fan = sqlite3.connect(fan_path)

        cell_lines = Database.load_tableList(con_fan)

        # output
        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_{}_{}.db'.format(ref, hbw))
        con_out = sqlite3.connect(out_path)

        M_ = len(cell_lines)

        tlist = Database.load_tableList(ref_con)
        for tname in tlist:
            chromosome, strand = tname.split('_')
            print(chromosome, strand)
            if chromosome == 'chrM' or chromosome == 'chrY':
                continue

            df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)
            if df_ref.empty:
                continue

            tss = df_ref[['start', 'end']].astype(int).mean(axis=1).astype(int)
            df_ref['start'] = tss.astype(int) - hbw
            df_ref['end'] = tss.astype(int) + hbw
            df_ref['corr'] = None

            contents = Parallel(n_jobs=num_cores)(delayed(self.processInput_corr)(df_ref, [chromosome, strand], cline, 100 * (i + 1) / M_)
                                                 for i, cline in enumerate(cell_lines))
            df_buf = pd.DataFrame(data=contents, index=cell_lines).T
            pd.concat([df_ref[['miRNA', 'start', 'end']], df_buf], axis=1).to_sql('_'.join([chromosome, strand]), con_out, index=None, if_exists='replace')

    def correlation_fan(self, hbw, ref='gene'):
        if ref == 'mir':
            # reference GENCODE
            ref_path = os.path.join(self.root, 'database', 'consistent_miRNA_330.db')
            ref_con = sqlite3.connect(ref_path)
        else:
            # reference GENCODE
            ref_path = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
            ref_con = sqlite3.connect(ref_path)

        # Fantom5
        fan_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_hCAGE_celllines.db')
        con_fan = sqlite3.connect(fan_path)

        cell_lines = Database.load_tableList(con_fan)

        # output
        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_{}_{}.db'.format(ref, hbw))
        con_out = sqlite3.connect(out_path)

        M_ = len(cell_lines)

        tlist = Database.load_tableList(ref_con)
        for tname in tlist:
            chromosome, strand = tname.split('_')
            print(chromosome, strand)
            if chromosome == 'chrM' or chromosome == 'chrY':
                continue

            df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)
            if df_ref.empty:
                continue

            if strand == '+':
                offset = 'start'
            else:
                offset = 'end'

            tss = deepcopy(df_ref[offset])
            df_ref['start'] = tss.astype(int) - hbw
            df_ref['end'] = tss.astype(int) + hbw
            df_ref['corr'] = None

            N_ = df_ref.shape[0]
            print(chromosome, strand)

            df_buf = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *cell_lines])
            df_buf[['start', 'end']] = df_ref[['start', 'end']]

            for i, cline in enumerate(cell_lines):
                df_fan = pd.read_sql_query("SELECT start, end, score FROM '{}' WHERE chromosome='{}' AND strand='{}'"
                                           "".format(cline, chromosome, strand), con_fan).sort_values(by=['start'])
                df_buf.loc[:, cline] = self.get_score_sum(df_ref[['start', 'end']], df_fan[['start', 'end', 'score']])
            pd.concat([df_ref['transcript_name'], df_buf], axis=1).to_sql('_'.join([chromosome, strand]), con_out, index=None, if_exists='replace')

    def get_rna_seq(self, df_ref, tname, con):
        df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
        df.index = df['start'].astype(str) + '-' + df['end'].astype(str)
        ref_index = df_ref['start'].astype(str) + '-' + df_ref['end'].astype(str)
        tns = set.intersection(set(ref_index), set(df.index))
        df_res = df.loc[tns]
        return df_res

    def correlation_fan_rna(self, hbw=500):
        # tissue list
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'tissues_fantom_rna.xlsx')
        df_tis = pd.read_excel(fpath, sheet_name='Sheet1')
        tissues = df_tis['RNA-seq']

        # reference GENCODE
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v32lift37.annotation_attr.db')
        ref_con = sqlite3.connect(ref_path)

        # Fantom5
        fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'FANTOM_tissue_spt.db')
        con_fan = sqlite3.connect(fan_path)

        # RNA-seq
        rna_path = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_spt.db')
        con_rna = sqlite3.connect(rna_path)

        # output
        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_{}.db'.format(hbw))
        con_out = sqlite3.connect(out_path)

        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_{}.db'.format(hbw))
        con_corr_out = sqlite3.connect(out_path)

        M_ = len(tissues)

        tlist = Database.load_tableList(ref_con)
        for tname in tlist:
            chromosome, strand = tname.split('_')
            if chromosome == 'chrM' or chromosome == 'chrY':
                continue

            sql = "SELECT start, end, gene_id, gene_name, gene_type, transcript_id, transcript_type, transcript_name " \
                  "FROM '{}' WHERE feature='transcript'"
            df_ref = pd.read_sql_query(sql.format(tname), ref_con)
            if df_ref.empty:
                continue

            if strand == '+':
                offset = 'start'
            else:
                offset = 'end'

            tss = deepcopy(df_ref[offset])
            df_ref['start'] = tss - hbw
            df_ref['end'] = tss + hbw
            df_ref['corr'] = None

            N_ = df_ref.shape[0]
            print(chromosome, strand)

            df_buf = {}
            for src in ['fantom', 'rna-seq']:
                df_buf[src] = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *tissues])
                df_buf[src][['start', 'end']] = df_ref[['start', 'end']]

            for i, tissue in enumerate(tissues):
                tname_rsc = '_'.join([tissue, chromosome, strand])
                df_fan = pd.read_sql_query("SELECT start, end, score FROM '{}' WHERE chromosome='{}' AND strand='{}'"
                                           "".format(tname_rsc, chromosome, strand), con_fan).sort_values(by=['start'])
                df_rna = pd.read_sql_query("SELECT start, end, FPKM FROM '{}' WHERE chromosome='{}' AND strand='{}'"
                                           "".format(tname_rsc, chromosome, strand), con_rna).sort_values(by=['start'])
                df_rna = df_rna[['start', 'end', 'FPKM']].rename(columns={"FPKM": "score"})

                for src, df_src in zip(['fantom', 'rna-seq'], [df_fan, df_rna]):
                    df_buf[src].loc[:, tissue] = self.get_score_sum(df_ref[['start', 'end']], df_src[['start', 'end', 'score']])

            for src in ['fantom', 'rna-seq']:
                df_buf[src].to_sql('_'.join([src, chromosome, strand]), con_out, index=None, if_exists='replace')

            for idx in df_ref.index:
                corr_coeff = np.corrcoef(df_buf['fantom'].loc[idx, tissues], df_buf['rna-seq'].loc[idx, tissues])
                df_ref.loc[idx, 'corr'] = corr_coeff[0, 1]
            df_ref.to_sql(tname, con_corr_out, index=None, if_exists='replace')

    def correlation(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna.db')
        con = sqlite3.connect(fpath)

        # reference GENCODE
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation_spt.db')
        ref_con = sqlite3.connect(ref_path)

        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'corr_fan_rna.db')
        con_out = sqlite3.connect(out_path)

        tlist = Database.load_tableList(ref_con)
        for tname in tlist:
            source, version, chromosome, strand = tname.split('.')
            df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)

            dfs_src = {}
            df_ref['corr'] = None
            for source in ['fantom', 'rna-seq']:
                tname_scr = '_'.join([source, chromosome, strand])
                dfs_src[source] = pd.read_sql_query("SELECT * FROM '{}'".format(tname_scr), con)

            for idx in df_ref.index:
                corr = np.corrcoef(dfs_src['fantom'].iloc[idx, 2:], dfs_src['rna-seq'].iloc[idx, 2:])
                df_ref.loc[idx, 'corr'] = corr[0, 1]
            df_ref.to_sql(tname, con_out, index=None, if_exists='replace')

    def high_correlation(self, hbw=500):
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        fpath = os.path.join(dirname, 'correlation_fan_rna_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        fpath = os.path.join(self.root, 'gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
        con_out = sqlite3.connect(fpath)
        dfs = []
        for tname in tlist:
            chromosome, strand = tname.split('_')
            df = pd.read_sql_query("SELECT * FROM '{}' WHERE corr>0.9 AND transcript_type='protein_coding'".format(tname), con)
            df.to_sql(tname, con_out, if_exists='replace', index=None)
            df.loc[:, 'chromosome'] = chromosome
            df.loc[:, 'strand'] = strand
            dfs.append(df)
        pd.concat(dfs).to_excel(fpath.replace('.db', '.xlsx'), index=None)


if __name__ == '__main__':
    cor = Correlation()
    if cor.hostname == 'mingyu-Precision-Tower-781':
        cor.to_server()
    else:
        from Regression import Regression

        rg = Regression()
        for hbw in [100, 300, 500]:
            cor.correlation_fan_rna(hbw)
            cor.high_correlation(hbw)
            cor.correlation_fan(hbw, ref='gene')
            cor.correlation_fan_cpu(hbw, ref='mir')
            rg.regression(hbw)
