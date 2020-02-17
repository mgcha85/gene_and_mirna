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

        server.job_script(fname, src_root=server_root, time='08:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def correlation_fan_rna_cpu(self, hbw=100, corr='spearman'):
        if corr == 'spearman':
            from scipy.stats import spearmanr
            corr = spearmanr
        else:
            from scipy.stats import pearsonr
            corr = pearsonr

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
                  "FROM '{}' WHERE feature='transcript' AND transcript_type='protein_coding'"
            df_ref = pd.read_sql_query(sql.format(tname), ref_con, index_col='transcript_id')
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
            df_ref['pvalue'] = None

            N_ = df_ref.shape[0]
            print(chromosome, strand)

            df_buf = {}
            for src in ['fantom', 'rna-seq']:
                df_buf[src] = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *tissues])
                df_buf[src][['start', 'end']] = df_ref[['start', 'end']]

            for i, tissue in enumerate(tissues):
                tname_rsc = '_'.join([tissue, chromosome, strand])
                df_fan = pd.read_sql_query("SELECT start, end, score FROM '{}'"
                                           "".format(tname_rsc, chromosome, strand), con_fan).sort_values(by=['start'])
                df_rna = pd.read_sql_query("SELECT start, end, FPKM, reference_id FROM '{}'"
                                           "".format(tname_rsc, chromosome, strand), con_rna, index_col='reference_id').sort_values(by=['start'])
                tid = set.intersection(set(df_ref.index), set(df_rna.index))
                df_buf['rna-seq'].loc[tid, tissue] = df_rna.loc[tid, 'FPKM']

                for src, df_src in zip(['fantom'], [df_fan]):
                    df_buf[src].loc[:, tissue] = self.get_score_sum(df_ref[['start', 'end']], df_src[['start', 'end', 'score']])

            df_buf = self.filtering(df_buf)
            for src in ['fantom', 'rna-seq']:
                df_buf[src].to_sql('_'.join([src, chromosome, strand]), con_out, if_exists='replace')

            for idx in df_buf['fantom'].index:
                rho, pvalue = corr(df_buf['fantom'].loc[idx, tissues], df_buf['rna-seq'].loc[idx, tissues])
                df_ref.loc[idx, 'corr'] = rho
                df_ref.loc[idx, 'pvalue'] = pvalue
            df_ref.dropna(subset=['corr']).to_sql(tname, con_corr_out, if_exists='replace')

    def correlation_fan_rna(self, hbw=0, corr='spearman'):
        if corr == 'spearman':
            from corr_gpu import Spearman
            corr = Spearman(self.root)
        else:
            from corr_gpu import Pearson
            corr = Pearson(self.root)

        # tissue list
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'tissues_fantom_rna.xlsx')
        df_tis = pd.read_excel(fpath, sheet_name='Sheet1')
        tissues = df_tis['RNA-seq']

        # reference CAGE clusters
        ref_path = os.path.join(self.root, 'database/Fantom/v5/cluster', 'result_all_clusters_spt.db')
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

            sql = "SELECT start, end, gene_name, transcript_id FROM '{}'"
            df_ref = pd.read_sql_query(sql.format(tname), ref_con, index_col='transcript_id')
            if df_ref.empty:
                continue
            df_ref['corr'] = None

            N_ = df_ref.shape[0]
            print(chromosome, strand)

            df_buf = {}
            for src in ['fantom', 'rna-seq']:
                df_buf[src] = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *tissues])
                df_buf[src][['start', 'end']] = df_ref[['start', 'end']]

            for i, tissue in enumerate(tissues):
                tname_rsc = '_'.join([tissue, chromosome, strand])
                df_fan = pd.read_sql_query("SELECT start, end, score FROM '{}'"
                                           "".format(tname_rsc, chromosome, strand), con_fan).sort_values(by=['start'])
                df_rna = pd.read_sql_query("SELECT start, end, FPKM, reference_id FROM '{}'"
                                           "".format(tname_rsc, chromosome, strand), con_rna, index_col='reference_id').sort_values(by=['start'])
                tid = set.intersection(set(df_ref.index), set(df_rna.index))
                df_buf['rna-seq'].loc[tid, tissue] = df_rna.loc[tid, 'FPKM'].astype(float)

                for src, df_src in zip(['fantom'], [df_fan]):
                    df_buf[src].loc[:, tissue] = self.get_score_sum(df_ref[['start', 'end']], df_src[['start', 'end', 'score']])

            df_buf = self.filtering(df_buf)
            for src in ['fantom', 'rna-seq']:
                df_buf[src].to_sql('_'.join([src, chromosome, strand]), con_out, if_exists='replace')

            df_ref['corr'] = pd.Series(data=corr.run(df_buf['fantom'][tissues], df_buf['rna-seq'][tissues], prod=False),
                                  index=df_buf['fantom'].index)
            df_ref.dropna(subset=['corr']).to_sql(tname, con_corr_out, if_exists='replace')

    def high_clusters(self, hbw):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_max.db'))
        tlist = Database.load_tableList(con)

        for tname in tlist:
            df = pd.read_sql("SELECT transcript_id, start, end, gene_name, "
                        "max(corr) FROM '{}' GROUP BY gene_name".format(tname), con)
            df = df.rename(columns={'max(corr)': 'corr'})
            df[['start', 'end']] = df[['start', 'end']].astype(int)
            df['corr'] = df['corr'].astype(float)
            df = df.sort_values(by=['start', 'end'])
            df.to_sql(tname, con_out, index=None, if_exists='replace')

    def high_correlation_by_thres(self, hbw=500):
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        fpath = os.path.join(dirname, 'correlation_fan_rna_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        contents = []
        xaxis = np.arange(0.5, 1, 0.1)
        for thres in xaxis:
            dfs = []
            for tname in tlist:
                chromosome, strand = tname.split('_')
                df = pd.read_sql_query("SELECT * FROM '{}' WHERE corr>{} AND transcript_type='protein_coding'".format(tname, thres), con)
                if df.empty:
                    continue

                df.loc[:, 'chromosome'] = chromosome
                df.loc[:, 'strand'] = strand
                dfs.append(df)
            df_res = pd.concat(dfs)
            contents.append(df_res.shape[0])
        plt.plot(xaxis, contents)
        plt.xlabel('correlation coefficient threshold')
        plt.ylabel('# protein coding')
        plt.grid()
        plt.show()

    def high_correlation(self, hbw=100, thres=0.8):
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        fpath = os.path.join(dirname, 'correlation_fan_rna_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        fpath = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
        con_out = sqlite3.connect(fpath)
        dfs = []
        cnt = 0
        for tname in tlist:
            chromosome, strand = tname.split('_')
            sql = "SELECT transcript_id, start, end, gene_name, MAX(corr) AS corr FROM (SELECT * FROM '{}' WHERE " \
                  "corr>{}) GROUP BY gene_name".format(tname, thres)
            # sql = "SELECT * FROM (SELECT transcript_id, start, end, gene_name, MAX(corr) AS corr FROM '{}' GROUP " \
            #       "BY gene_name) WHERE corr>{}".format(tname, thres)
            df = pd.read_sql_query(sql, con)
            if df.empty:
                continue

            df.to_sql(tname, con_out, if_exists='replace', index=None)
            df.loc[:, 'chromosome'] = chromosome
            df.loc[:, 'strand'] = strand
            cnt += df.shape[0]
            dfs.append(df)

        if len(dfs) > 0:
            print(cnt)
            pd.concat(dfs).to_excel(fpath.replace('.db', '.xlsx'), index=None)

    def sum_fan(self, hbw, ref='gene'):
        if ref == 'mir':
            # reference consistent miRNA
            ref_path = os.path.join(self.root, 'database', 'consistent_miRNA_330_spt.db')
            ref_con = sqlite3.connect(ref_path)
            label = 'miRNA'
        else:
            # reference CAGE clusters
            # ref_path = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
            ref_path = os.path.join(self.root, 'database/gencode', 'other_researches_inter_fan_rna_{}.db'.format(hbw))
            ref_con = sqlite3.connect(ref_path)
            label = 'transcript_id'

        # Fantom5
        fan_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_hCAGE_celllines.db')
        # fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'FANTOM_tissue.db')
        con_fan = sqlite3.connect(fan_path)

        cell_lines = Database.load_tableList(con_fan)

        # output
        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_{}_{}_others_inter.db'.format(ref, hbw))
        # out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_{}_{}.db'.format(ref, hbw))
        con_out = sqlite3.connect(out_path)

        M_ = len(cell_lines)
        for tname in Database.load_tableList(ref_con):
            chromosome, strand = tname.split('_')
            print(chromosome, strand)
            if chromosome == 'chrM' or chromosome == 'chrY':
                continue

            df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)
            if df_ref.empty:
                continue

            df_ref['corr'] = None
            if ref == 'mir':
                df_ref['start'] = df_ref['tss'] - hbw
                df_ref['end'] = df_ref['tss'] + hbw

            N_ = df_ref.shape[0]
            df_buf = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *cell_lines])
            df_buf[['start', 'end']] = df_ref[['start', 'end']]

            for i, cline in enumerate(cell_lines):
                df_fan = pd.read_sql_query("SELECT start, end, score FROM '{}' WHERE chromosome='{}' AND strand='{}'"
                                           "".format(cline, chromosome, strand), con_fan).sort_values(by=['start'])
                print(df_ref.shape[0], df_fan.shape[0])
                if df_fan.shape[0] == 0:
                    continue
                df_buf.loc[:, cline] = self.get_score_sum(df_ref[['start', 'end']], df_fan[['start', 'end', 'score']])
            pd.concat([df_ref[label], df_buf], axis=1).to_sql('_'.join([chromosome, strand]), con_out, index=None, if_exists='replace')

    # def sum_fan(self, hbw, ref='gene'):
    #     if ref == 'mir':
    #         # reference consistent miRNA
    #         ref_path = os.path.join(self.root, 'database', 'consistent_miRNA_330_spt.db')
    #         ref_con = sqlite3.connect(ref_path)
    #     else:
    #         # reference GENCODE
    #         ref_path = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
    #         ref_con = sqlite3.connect(ref_path)
    # 
    #     # Fantom5
    #     fan_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_hCAGE_celllines.db')
    #     con_fan = sqlite3.connect(fan_path)
    # 
    #     cell_lines = Database.load_tableList(con_fan)
    # 
    #     # output
    #     out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_{}_{}.db'.format(ref, hbw))
    #     con_out = sqlite3.connect(out_path)
    # 
    #     M_ = len(cell_lines)
    # 
    #     tlist = Database.load_tableList(ref_con)
    #     for tname in tlist:
    #         chromosome, strand = tname.split('_')
    #         print(chromosome, strand)
    #         if chromosome == 'chrM' or chromosome == 'chrY':
    #             continue
    # 
    #         df_ref = pd.read_sql_query("SELECT * FROM '{}'".format(tname), ref_con)
    #         if df_ref.empty:
    #             continue
    # 
    #         if strand == '+':
    #             offset = 'start'
    #         else:
    #             offset = 'end'
    # 
    #         tss = deepcopy(df_ref[offset])
    #         df_ref['start'] = tss.astype(int) - hbw
    #         df_ref['end'] = tss.astype(int) + hbw
    #         df_ref['corr'] = None
    # 
    #         N_ = df_ref.shape[0]
    #         print(chromosome, strand)
    # 
    #         df_buf = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *cell_lines])
    #         df_buf[['start', 'end']] = df_ref[['start', 'end']]
    # 
    #         for i, cline in enumerate(cell_lines):
    #             df_fan = pd.read_sql_query("SELECT start, end, score FROM '{}' WHERE chromosome='{}' AND strand='{}'"
    #                                        "".format(cline, chromosome, strand), con_fan).sort_values(by=['start'])
    #             df_buf.loc[:, cline] = self.get_score_sum(df_ref[['start', 'end']], df_fan[['start', 'end', 'score']])
    #         pd.concat([df_ref['transcript_name'], df_buf], axis=1).to_sql('_'.join([chromosome, strand]), con_out, index=None, if_exists='replace')

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

    def get_rna_seq(self, df_ref, tname, con):
        df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
        df.index = df['start'].astype(str) + '-' + df['end'].astype(str)
        ref_index = df_ref['start'].astype(str) + '-' + df_ref['end'].astype(str)
        tns = set.intersection(set(ref_index), set(df.index))
        df_res = df.loc[tns]
        return df_res

    def filtering(self, dfs):
        ridx = []
        for idx in dfs['fantom'].index:
            fan_row = dfs['fantom'].loc[idx, 'appendix':].astype(float)
            rna_row = dfs['rna-seq'].loc[idx, 'appendix':].astype(float)
            fmidx = fan_row.idxmax()
            rmidx = rna_row.idxmax()

            if fan_row[fmidx] == 0 or rna_row[rmidx] == 0:
                ridx.append(idx)
            elif fan_row.median() == 0 or rna_row.median() == 0:
                ridx.append(idx)
            elif fan_row[fmidx] > fan_row.drop(fmidx).sum() or rna_row[fmidx] > rna_row.drop(fmidx).sum():
                ridx.append(idx)

        for type in ['fantom', 'rna-seq']:
            dfs[type] = dfs[type].drop(ridx)
        return dfs

    def correlation_gpu(self, hbw, opt, corr='spearman'):
        if corr == 'spearman':
            from corr_gpu import Spearman
            corr = Spearman(self.root)
        else:
            from corr_gpu import Pearson
            corr = Pearson(self.root)

        fpath_mir = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_mir_100.db')
        con_mir = sqlite3.connect(fpath_mir)

        fpath_gene = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_gene_{}_others_inter.db'.format(hbw))
        con_gene = sqlite3.connect(fpath_gene)

        fpath_out = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}_{}_others.db'.format(hbw, opt))
        con_out = sqlite3.connect(fpath_out)

        def merge(con):
            tlist = Database.load_tableList(con)
            dfs = []
            for tname in tlist:
                dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
            return pd.concat(dfs).set_index(dfs[-1].columns[0])

        df_mir = merge(con_mir)
        df_gene = merge(con_gene)
        df_out = pd.DataFrame(data=corr.run(df_gene.loc[:, '10964C':], df_mir.loc[:, '10964C':]), index=df_gene.index, columns=df_mir.index)
        df_out.to_sql('corr', con_out, if_exists='replace')

    def get_sample_corr(self, hbw):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)

        df = pd.read_excel(os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.xlsx'.format(hbw)))
        ridx = np.random.randint(0, df.shape[0], 3)
        df = df.iloc[ridx, :]

        contents = []
        for idx in df.index:
            chromosome = df.loc[idx, 'chromosome']
            strand = df.loc[idx, 'strand']
            start = df.loc[idx, 'start']
            end = df.loc[idx, 'end']

            df_fan = pd.read_sql("SELECT * FROM 'fantom_{}_{}' WHERE start={} AND end={}".format(chromosome, strand, start, end), con)
            df_rna = pd.read_sql("SELECT * FROM 'rna-seq_{}_{}' WHERE start={} AND end={}".format(chromosome, strand, start, end), con)
            fan = df_fan.loc[0, "appendix":].astype(float)
            rna = df_rna.loc[0, "appendix":].astype(float)
            contents.append([chromosome, start, end, strand, ','.join(fan.round(4).astype(str)), ','.join(rna.round(4).astype(str))])
        df_res = pd.DataFrame(data=contents, columns=['chromosome', 'start', 'end', 'strand', 'FANTOM', 'RNA-seq'])
        df_res.to_excel(os.path.join(self.root, 'database/Fantom/v5/tissues', 'corr_sample.xlsx'), index=None)

    def plot_sample(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'corr_sample.xlsx')
        df = pd.read_excel(fpath)

        for i, idx in enumerate(df.index):
            fantom = np.array(list(map(float, df.loc[idx, 'FANTOM'].split(','))))
            rna_seq = np.array(list(map(float, df.loc[idx, 'RNA-seq'].split(','))))
            plt.subplot(3, 2, 2*(i+1)-1)
            plt.scatter(np.arange(len(fantom)), fantom)
            plt.title('FANTOM (sample{})'.format(i+1))
            plt.xlabel('tissues')
            plt.ylabel('expression value')
            plt.grid()

            plt.subplot(3, 2, 2*(i+1))
            plt.scatter(np.arange(len(fantom)), rna_seq)
            plt.title('RNA-seq (sample{})'.format(i+1))
            plt.xlabel('tissues')
            plt.ylabel('expression value')
            plt.grid()
        plt.show()

    def add_corr(self, hbw):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'result'", con, index_col='miRNA')
        transcripts = df['Transcripts'].str.split(';')

        for i, mir in enumerate(df.index):
            print('{:,d} / {:,d}'.format(i+1, df.shape[0]))
            tname = transcripts[mir]

            sql = ["transcript_id='{}'".format(x) for x in tname]
            sql = ' OR '.join(sql)
            sql = "SELECT * FROM 'corr' WHERE " + sql
            corr_coef = pd.read_sql(sql, con, index_col='transcript_id')
            if corr_coef.empty:
                continue
            df.loc[mir, 'corr'] = ';'.join(corr_coef[mir].round(4).astype(str).values)
            stats = np.array([corr_coef[mir].mean(), corr_coef[mir].std(), np.median(corr_coef[mir]),
                      corr_coef[mir].min(), corr_coef[mir].max()])
            df.loc[mir, 'corr_stats'] = ';'.join(stats.round(4).astype(str))
            df.loc[mir, 'n'] = corr_coef[mir].shape[0]
        df.to_sql('result', con, if_exists='replace')


if __name__ == '__main__':
    cor = Correlation()
    if cor.hostname == 'mingyu-Precision-Tower-7810':
        cor.to_server()
    else:
        from Regression import Regression
        from mir_gene import mir_gene
        from set_go import set_go

        rg = Regression()
        mg = mir_gene()
        sg = set_go()

        # cor.high_correlation_by_thres(100)
        # cor.get_sample_corr(hbw)
        # cor.plot_sample()

        for hbw in [100]:
            for opt in ['nz']:
                # cor.correlation_fan_rna(hbw)
                # cor.high_clusters(hbw)
                # cor.high_correlation(hbw, 0.8)

                # cell lines
                cor.sum_fan(hbw, ref='gene')
                # cor.sum_fan(hbw, ref='mir')

                # rg.regression(hbw, opt)
                # rg.report(hbw, opt)
                # rg.add_gene_name(hbw, opt)
                # rg.filtering(hbw)

                cor.correlation_gpu(hbw, opt)
                # cor.add_corr(hbw)

                # mg.comparison(hbw)
                # mg.phypher(hbw)
                # mg.plot(hbw)

                # sg.set_input(hbw, opt)
                # sg.submit_data(hbw, opt, bg=True)
                # sg.extract_genes(hbw, opt)

                # sg.result(hbw, opt)
                # sg.to_tg(hbw, opt)

                # sg.move(hbw, opt)
                # sg.hyper_test(hbw, opt)
                # sg.plot_hyper_test(hbw, opt)
