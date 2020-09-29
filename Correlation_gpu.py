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
from scipy.stats import spearmanr


class Correlation:
    def __init__(self, root):
        self.root = root

    def to_server(self, root, tag):
        from Server import Server
        import sys

        which = 'newton'
        server = Server(root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname__ = os.path.split(local_path)
        fname = fname__.replace('.py', '_{}.py'.format(tag))

        curdir = os.getcwd().split(os.sep)[-1]
        server_root = os.path.join(server.server, 'source', curdir).replace(os.sep, '/')
        server_path = local_path.replace(dirname, server_root)
        server_path = server_path.replace(fname__, fname)

        server.job_script(fname, src_root=server_root, time='08:00:00', pversion=3)
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', server_root + '/' + 'dl-submit.slurm')

        stdin, stdout, stderr = server.ssh.exec_command("cd {root};sbatch {root}/dl-submit.slurm".format(root=server_root))
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
            df_ref = pd.read_sql(sql.format(tname), ref_con, index_col='transcript_id')
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
                df_fan = pd.read_sql("SELECT start, end, score FROM '{}'"
                                           "".format(tname_rsc, chromosome, strand), con_fan).sort_values(by=['start'])
                df_rna = pd.read_sql("SELECT start, end, FPKM, reference_id FROM '{}'"
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
        # ref_path = os.path.join(self.root, 'database/Fantom/v5/cluster', 'result_all_clusters_spt.db')
        # ref_con = sqlite3.connect(ref_path)
        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v32lift37.annotation_attr_spt.db')
        ref_con = sqlite3.connect(ref_path)

        # # for check (temp)
        # with open('remain.txt', 'rt') as f:
        #     remains = f.read().split('\n')

        # Fantom5hg19.cage_peak_phase1and2combined_tpm_ann.osc_avg
        fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'CAGE_tag_tissue_spt.db')
        # fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'FANTOM_tissue_spt.db')
        con_fan = sqlite3.connect(fan_path)

        # RNA-seq
        rna_path = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_spt.db')
        con_rna = sqlite3.connect(rna_path)

        # output
        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_{}_processed.db'.format(hbw))
        con_out = sqlite3.connect(out_path)

        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_{}_processed.db'.format(hbw))
        con_corr_out = sqlite3.connect(out_path)

        M_ = len(tissues)

        dfs = []
        tlist = Database.load_tableList(ref_con)
        for tname in tlist:
            chromosome, strand = tname.split('_')
            # if chromosome == 'chrM' or chromosome == 'chrY':
            #     continue

            sql = "SELECT start, end, gene_name, transcript_id FROM '{}'"
            df_ref = pd.read_sql(sql.format(tname), ref_con, index_col='transcript_id')
            if df_ref.empty:
                continue

            # # for check (temp)
            # df_ref = df_ref.loc[set.intersection(set(df_ref.index), set(remains))]
            # if df_ref.empty:
            #     continue

            df_ref['corr'] = None

            N_ = df_ref.shape[0]
            print(chromosome, strand)

            df_buf = {}
            for src in ['fantom', 'rna-seq']:
                df_buf[src] = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *tissues])
                df_buf[src][['start', 'end']] = df_ref[['start', 'end']]

            for i, tissue in enumerate(tissues):
                tname_rsc = '_'.join([tissue, chromosome, strand])
                # RNA-seq
                df_rna = pd.read_sql("SELECT start, end, FPKM, reference_id FROM '{}'"
                                           "".format(tname_rsc, chromosome, strand), con_rna, index_col='reference_id').sort_values(by=['start'])
                df_rna['FPKM'] = df_rna['FPKM'].round(4)
                tid = set.intersection(set(df_ref.index), set(df_rna.index))
                df_buf['rna-seq'].loc[tid, tissue] = df_rna.loc[tid, 'FPKM'].astype(float)

                # FANTOM
                df_fan = pd.read_sql("SELECT start, end, score FROM '{}'"
                                           "".format(tname_rsc, chromosome, strand), con_fan).sort_values(by=['start'])
                df_fan['score'] = df_fan['score'].round(4)

                for src, df_src in zip(['fantom'], [df_fan]):
                    df_buf[src].loc[:, tissue] = self.get_score_sum(df_ref[['start', 'end']], df_src[['start', 'end', 'score']])

            for src in ['fantom', 'rna-seq']:
                df_buf[src].to_sql('_'.join([src, chromosome, strand, '_raw']), con_out, if_exists='replace')
            df_buf = self.filtering(df_buf)
            for src in ['fantom', 'rna-seq']:
                df_buf[src].to_sql('_'.join([src, chromosome, strand]), con_out, if_exists='replace')

            df_ref['corr'] = pd.Series(data=corr.run(df_buf['fantom'][tissues].round(4), df_buf['rna-seq'][tissues].round(4), prod=False),
                                  index=df_buf['fantom'].index).round(2)
            df_ref = df_ref.dropna(subset=['corr'])
            df_ref.to_sql(tname, con_corr_out, if_exists='replace')
            dfs.append(df_ref)
        pd.concat(dfs).to_excel(out_path.replace('.db', '.xlsx'))

    def corr_stats(self, hbw):
        import pandas_profiling
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_{}.xlsx'.format(hbw))
        df = pd.read_excel(fpath, index_col=0)

        mean = df['corr'].mean()
        median = df['corr'].median()
        max = df['corr'].max()
        min = df['corr'].min()
        std = df['corr'].std()
        quantile7 = df['corr'].quantile(q=0.7)
        quantile8 = df['corr'].quantile(q=0.8)
        quantile9 = df['corr'].quantile(q=0.9)
        print('mean: {:0.2f}, median: {:0.2f}, max: {:0.2f}, min: {:0.2f}, std: {:0.2f}, quantile7: {:0.2f}, '
              'quantile8: {:0.2f}, quantile9: {:0.2f}'
              ''.format(mean, median, max, min, std, quantile7, quantile8, quantile9))
        
    def high_clusters(self, hbw):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_{}_processed.db'.format(hbw))
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
                df = pd.read_sql("SELECT * FROM '{}' WHERE corr>{} AND transcript_type='protein_coding'".format(tname, thres), con)
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
        fpath = os.path.join(dirname, 'correlation_fan_rna_{}_processed.db'.format(hbw))
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        fpath = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
        con_out = sqlite3.connect(fpath)
        dfs = []
        cnt = 0
        for tname in tlist:
            chromosome, strand = tname.split('_')
            print(chromosome, strand)
            # sql = "SELECT transcript_id, start, end, gene_name, MAX(corr) AS corr FROM (SELECT * FROM '{}' WHERE " \
            #       "round(corr, 2)>={}) GROUP BY gene_name".format(tname, thres)
            # sql = "SELECT * FROM '{}' WHERE corr>{}".format(tname, thres)

            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            if df.empty:
                continue

            df_grps = []
            for gname, df_grp in df.groupby('gene_name'):
                idx = df_grp['corr'].idxmax()
                df_grps.append(df_grp.loc[idx])
            df = pd.concat(df_grps, axis=1).T
            df = df[df['corr'] >= thres]
            if df.empty:
                continue

            df[['start', 'end']] = df[['start', 'end']].astype(int)
            df.to_sql(tname, con_out, if_exists='replace', index=False)
            df.loc[:, 'chromosome'] = chromosome
            df.loc[:, 'strand'] = strand
            cnt += df.shape[0]
            dfs.append(df)

        if len(dfs) > 0:
            print(cnt)
            pd.concat(dfs).to_excel(fpath.replace('.db', '.xlsx'), index=None)

    def sum_fan(self, hbw, ref='gene', type='cell_lines'):
        if ref == 'mir':
            # reference consistent miRNA
            ref_path = os.path.join(self.root, 'database', 'consistent_miRNA_330_spt.db')
            ref_con = sqlite3.connect(ref_path)
            label = 'miRNA'
        else:
            # reference CAGE clusters
            ref_path = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
            # ref_path = os.path.join(self.root, 'database/gencode', 'other_researches_union_fan_rna_{}.db'.format(hbw))
            ref_con = sqlite3.connect(ref_path)
            label = 'transcript_id'

        # Fantom5
        if type == 'cell_lines':
            fan_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'hg19.cage_peak_phase1and2combined_tpm_ann.osc_avg.db')
        elif type == 'tissues':
            fan_path = os.path.join(self.root, 'database/Fantom/v5/tissues', 'CAGE_tag_tissue_spt.db')
        else:
            raise('wrong type')

        con_fan = sqlite3.connect(fan_path)
        out_path = os.path.join(self.root, 'database/Fantom/v5/{}'.format(type), 'sum_fan_{}_{}.db'.format(ref, hbw))

        cell_lines = sorted(list(set([x.split('_')[0] for x in Database.load_tableList(con_fan)])))
        # output
        con_out = sqlite3.connect(out_path)

        M_ = len(cell_lines)
        for tname in Database.load_tableList(ref_con):
            chromosome, strand = tname.split('_')
            print(chromosome, strand)
            if chromosome == 'chrM' or chromosome == 'chrY':
                print('[SKIP] chromosome is {}'.format(chromosome))
                continue

            df_ref = pd.read_sql("SELECT * FROM '{}'".format(tname), ref_con)
            if df_ref.empty:
                print('[SKIP] reference is empty')
                continue

            df_ref['corr'] = None
            if ref == 'mir':
                df_ref['start'] = df_ref['tss'] - hbw
                df_ref['end'] = df_ref['tss'] + hbw

            N_ = df_ref.shape[0]
            df_buf = pd.DataFrame(data=np.zeros((N_, M_ + 2)), index=df_ref.index, columns=['start', 'end', *cell_lines])
            df_buf[['start', 'end']] = df_ref[['start', 'end']]

            for i, cline in enumerate(cell_lines):
                # print(cline)
                if type == 'cell_lines':
                    df_fan = pd.read_sql("SELECT start, end, avg_score as score FROM '{}' WHERE chromosome='{}' AND strand='{}'"
                                               "".format(cline, chromosome, strand), con_fan).sort_values(by=['start'])
                else:
                    df_fan = pd.read_sql("SELECT start, end, score FROM '{}_{}_{}'"
                                               "".format(cline, chromosome, strand), con_fan).sort_values(by=['start'])
                if df_fan.shape[0] == 0:
                    print('[SKIP] fantom data is empty')
                    continue
                df_buf.loc[:, cline] = self.get_score_sum(df_ref[['start', 'end']], df_fan[['start', 'end', 'score']])
            pd.concat([df_ref[label], df_buf], axis=1).to_sql('_'.join([chromosome, strand]), con_out, index=False, if_exists='replace')

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
    #         df_ref = pd.read_sql("SELECT * FROM '{}'".format(tname), ref_con)
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
    #             df_fan = pd.read_sql("SELECT start, end, score FROM '{}' WHERE chromosome='{}' AND strand='{}'"
    #                                        "".format(cline, chromosome, strand), con_fan).sort_values(by=['start'])
    #             df_buf.loc[:, cline] = self.get_score_sum(df_ref[['start', 'end']], df_fan[['start', 'end', 'score']])
    #         pd.concat([df_ref['transcript_name'], df_buf], axis=1).to_sql('_'.join([chromosome, strand]), con_out, index=None, if_exists='replace')

    def get_score_sum(self, df_ref, df_res):
        from sqlgpu import Sqlgpu

        sqlgpu = Sqlgpu()
        out = sqlgpu.bin_run(df_ref, df_res)
        return out

    def get_score_sum_cpu(self, df_ref, df_res):
        df_out = pd.Series(index=df_ref.index)
        for idx in df_ref.index:
            if idx != 'ENST00000525827.5_1':
                continue
            start = df_ref.loc[idx, 'start']
            end = df_ref.loc[idx, 'end']
            df_out[idx] = df_res.loc[df_res[(df_res['start'] < end) & (df_res['end'] > start)].index, 'score'].sum()
        return df_out

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
            df = pd.read_sql(sql, con_fan).sort_values(by=['start'])
            contents[i] = df['score'].sum()
        return contents

    def get_rna_seq(self, df_ref, tname, con):
        df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
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

        fpath_gene = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'sum_fan_gene_{}.db'.format(hbw))
        con_gene = sqlite3.connect(fpath_gene)

        fpath_out = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}_{}.db'.format(hbw, opt))
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
        df_out = df_out.replace(2, np.nan)

        # df_out = pd.DataFrame(index=df_gene.index, columns=df_mir.index)
        # df_gene_tr = df_gene.loc[:, '10964C':].T
        # df_mir_tr = df_mir.loc[:, '10964C':].T
        # for mir in df_mir_tr.columns:
        #     df_corr = df_mir_tr[mir:mir].corrwith(df_gene_tr, method='spearman')
        #     print('see')
        # df_corr.to_excel('temp.xlsx')
        # exit(1)
        #
        # for i in range(df_corr.shape[0]):
        #     for j in range(df_corr.shape[1]):
        #         if i != j:
        #             gene = df_corr.index[i]
        #             mir = df_corr.columns[j]
        #             df_out.loc[gene, mir] = df_corr.iloc[i, j]

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
    hostname = socket.gethostname()
    if hostname == 'mingyu-Precision-Tower-7810':
        root = '/home/mingyu/Bioinformatics'
    elif hostname == 'DESKTOP-DLOOJR6' or hostname == 'DESKTOP-1NLOLK4':
        root = 'D:/Bioinformatics'
    elif hostname == 'mingyu-Inspiron-7559':
        root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
    else:
        root = '/lustre/fs0/home/mcha/Bioinformatics'

    cor = Correlation(root)
    if hostname == 'DESKTOP-DLOOJR6' or hostname == 'DESKTOP-1NLOLK4':
        cor.to_server(root, "")
    else:
        from Regression import Regression
        from mir_gene import mir_gene
        from set_go import set_go
        from validation import validation

        rg = Regression(root)
        mg = mir_gene()
        sg = set_go(root)
        val = validation(root)

        # cor.high_correlation_by_thres(100)
        # cor.get_sample_corr(hbw)
        # cor.plot_sample()

        for hbw in [100]:
            for opt in ['nz']:
                # cor.corr_stats(hbw)
                # cor.correlation_fan_rna(hbw)
                # cor.high_clusters(hbw)
                # cor.high_correlation(hbw, 0.75)
                # exit(1)

                for type in ['tissues']:
                    cor.sum_fan(hbw, ref='gene', type=type)
                    cor.sum_fan(hbw, ref='mir', type=type)
                    # rg.regression(hbw, opt, type)
                    # rg.report(hbw, opt, type)
                    rg.add_gene_name(hbw, opt, type)
                # rg.filtering(hbw)

                # cor.correlation_gpu(hbw, opt)
                # cor.add_corr(hbw)

                # mg.comparison(hbw)
                # mg.phypher(hbw)
                # mg.plot(hbw)

                # sg.set_input(hbw, opt)
                # sg.submit_data(hbw, opt, bg=True)
                # sg.extract_genes(hbw, opt)
                #
                # sg.result(hbw, opt)
                # sg.to_tg(hbw, opt)

                # sg.move(hbw, opt)
                # sg.hyper_test(hbw, opt)
                # sg.plot_hyper_test(hbw, opt)

                # rg.cross_regression(hbw, opt)
                # rg.cross_stats(100, opt)

                # val.wilcox(opt)
