import pandas as pd
import numpy as np
import sqlite3
import os
import socket
import sys

from Database import Database
from Server import Server
from joblib import Parallel, delayed
import multiprocessing
from scipy.stats import pearsonr
from joblib import Parallel, delayed
import multiprocessing
import pickle as pkl


class Correlation2:
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
        self.cells = []
        self.version = 0
        self.band = 100
        self.csize = 100
        self.num_cores = multiprocessing.cpu_count()

    def to_server(self):
        which = 'newton'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='12:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def get_reference(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql_query("SELECT * FROM 'human_genes'", con)
        df_unique = df.drop_duplicates(subset=['gene_name'], keep=False)

        idx_duplicates = list(set(df.index) - set(df_unique.index))
        df_duplicates = df.loc[idx_duplicates]
        df_duplicates = df_duplicates[df_duplicates['source'] == 'ENSEMBL']
        df_duplicates = df_duplicates.drop_duplicates(subset=['gene_name'], keep=False)
        return pd.concat([df_unique, df_duplicates]).reset_index(drop=True)

    def pre_data(self):
        fpath = os.path.join(self.root, 'software/stringtie-1.3.6', 'genes.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        out_con = sqlite3.connect(fpath.replace('.db', '_2.db'))

        for tname in tlist:
            df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            attribute = df_rna['attribute'].str.split('; ')

            pkm = []
            for attr in attribute:
                dict = {}
                for a in attr[:-1]:
                    key, value = a.split(' ')
                    dict[key] = value.replace('"', '')

                if 'gene_name' not in dict:
                    pkm.append(np.nan)
                    continue
                pkm.append(dict['gene_name'])

            df_res = pd.DataFrame(data=pkm, columns=['gene_name'])
            df_rna = pd.concat([df_rna, df_res], axis=1)
            df_rna.drop(['attribute', 'frame'], axis=1).to_sql(tname, out_con, index=None, if_exists='append')

    def processInput(self, df, gene, threshold, mean, std):
        chromosome = df['chromosome'].iloc[0]
        strand = df['strand'].iloc[0]

        if self.version == 1:
            df_sub = self.merge_tags(df, distance=self.csize)
            tss = ';'.join(df_sub['region'].astype(str))
            scores = ';'.join(df_sub['score'].astype(str))
            widths = ';'.join(df_sub['width'].astype(str))
            score = df_sub['score'].mean()
            score_norm = ';'.join(((df_sub['score'] - mean) / std).round(4).astype(str))

        elif self.version == 2:
            df_sub = self.merge_tags(df, distance=self.csize, threshold=threshold)
            tss = ';'.join(df_sub['region'].astype(str))
            scores = ';'.join(df_sub['score'].astype(str))
            widths = ';'.join(df_sub['width'].astype(str))
            score = df_sub['score'].mean()
            score_norm = ';'.join(((df_sub['score'] - mean) / std).round(4).astype(str))

        elif self.version == 3:
            midx = df['score'].idxmax()
            tss = str(df.loc[midx, 'start'])
            scores = str(df.loc[midx, 'score'])
            widths = str(df.loc[midx, 'end'] - df.loc[midx, 'start'])
            score_norm = str(df.loc[midx, 'score (norm)'])
            score = scores

        else:
            tss = ';'.join(df['start'].astype(str))
            scores = ';'.join(df['score'].astype(str))
            widths = ';'.join((df['end'] - df['start']).astype(str))
            score = df['score'].mean()
            score_norm = ';'.join(df.loc[:, 'score (norm)'].round(4).astype(str))

        return [chromosome, strand, tss, widths, scores, score, score_norm, gene]

    def fantom_unique_gene(self):
        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}.db'.format(self.band))
        con_fan = sqlite3.connect(fpath_fan)
        con_out = sqlite3.connect(fpath_fan.replace('.db', '_{}_v{}.db'.format(self.csize, self.version)))
        tlist_rna = Database.load_tableList(con_fan)

        if self.num_cores > 6:
            self.num_cores = 6
        print('num_cores: {}'.format(self.num_cores))

        threshold = None
        for tname in tlist_rna:
            print(tname)
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_fan)
            mean = df_fan['score'].mean()
            std = df_fan['score'].std()
            df_fan.loc[:, 'score (norm)'] = ((df_fan['score'] - mean) / std).round(4)
            if self.version == 2:
                threshold = self.get_threshold(df_fan)
            df_grp = df_fan.groupby('gene_name')

            contents = Parallel(n_jobs=self.num_cores)(delayed(self.processInput)(df_sub, gene, threshold, mean, std) for gene, df_sub in df_grp)
            df_res = pd.DataFrame(contents, columns=['chromosome', 'strand', 'tss', 'widths', 'scores', 'score', 'score_norm', 'gene_name'])
            df_res.to_sql(tname, con_out, if_exists='replace', index=None)

    def rna_unique_gene(self):
        fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_out.db')
        con_rna = sqlite3.connect(fpath_rna)
        con_out = sqlite3.connect(fpath_rna.replace('.db', '2.db'))
        tlist_rna = Database.load_tableList(con_rna)

        for tname in tlist_rna:
            print(tname)
            df_rna = pd.read_sql_query("SELECT * FROM {} WHERE FPKM>0".format(tname), con_rna)
            df_rna['FPKM'] = df_rna['FPKM'].astype(float)
            df_grp = df_rna.groupby('gene_name')

            contents = []
            for gene, df_sub in df_grp:
                midx = df_sub['FPKM'].idxmax()
                contents.append(df_sub.loc[midx:midx, :])

            df_res = pd.concat(contents)
            mean = df_res['FPKM'].mean()
            std = df_res['FPKM'].std()
            df_res.loc[:, 'FPKM (norm)'] = ((df_res['FPKM'] - mean) / std).round(4)
            df_res.to_sql(tname, con_out, if_exists='replace', index=None)

    def get_valid_tissues(self, tlist_rna, tlist_fan):
        table_names = {}
        for tname_fan in tlist_fan:
            tlist_rna_sub = [x for x in tlist_rna if tname_fan in x]
            table_names[tname_fan] = tlist_rna_sub
        return table_names

    def correlation(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'correlation_{}_{}_v{}.db'.format(self.band, self.csize, self.version))
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        dfs = {}
        gene_names = []
        for tname in tlist:
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con, index_col='gene_name')
            if df.empty:
                continue
            dfs[tname] = df
            gene_names.append(set(df.index))

        report = []
        gene_names = sorted(list(set.intersection(*gene_names)))
        for i, gname in enumerate(gene_names):
            if i % 100 == 0:
                print('{:0.2f}%'.format(100 * (i + 1) / len(gene_names)))

            scores = []
            fpkms = []
            tissues = []
            tsss = []
            widths = []
            num_cls = 0
            loci = ''
            tissue = ''
            for j, (tname, df) in enumerate(dfs.items()):
                tissue__ = tname.split('_')[0]
                if tissue != tissue__:
                    tissue = tissue__
                    tissues.append(tissue)
                    scores.append(df.loc[gname, 'score'])
                    fpkms.append(df.loc[gname, 'FPKM'])
                    tsss.append(df.loc[gname, 'tss'])
                    width = df.loc[gname, 'widths'].split(';')
                    num_cls += len(width)
                    widths.append(df.loc[gname, 'widths'])
                if j == 0:
                    chromosome = df.loc[gname, 'chromosome']
                    strand = df.loc[gname, 'strand']
                    if strand == '+':
                        tss = df.loc[gname, 'start']
                    else:
                        tss = df.loc[gname, 'end']
                    loc = ':'.join([chromosome, str(tss), strand])

            if i == 0:
                print(tissues)

            fpkms = np.array(fpkms)
            fpkm_str = fpkms.astype(str)
            scores = np.array(scores).astype(float)
            # num_cls = len(scores)
            num_cls /= len(scores)
            spearman_xcor = pearsonr(scores, fpkms)
            spec = '{:0.2f}:{}:{:0.4f}'.format(num_cls, ','.join(widths), scores.mean())

            # if np.isnan(spearman_xcor.correlation):
            #     continue
            report.append([gname, loc, ':'.join(tsss), spec,
                           ':'.join(scores.astype(str)), ','.join(fpkm_str), '{:0.4f}'.format(spearman_xcor[0]),
                           '{:0.4f}'.format(spearman_xcor[1])])

        df_rep = pd.DataFrame(data=report, columns=['gnames', 'loc', 'tss', 'spec', 'scores', 'fpkms', 'corr', 'p-value'])
        df_rep.to_excel('correlation_{}_{}_v{}.xlsx'.format(self.band, self.csize, self.version), index=None)

    def cut_by_threshold(self, df):
        for idx in df.index:
            start, end = df.loc[idx, 'region'][1:-1].split(',')

            score = df['score'].sort_values().drop_duplicates(keep='first')
            thres = pd.qcut(score, 10).iloc[0].right
            print(thres)
        return df[df['score'] > thres]

    def get_threshold(self, df):
        score = df['score'].sort_values()
        # score = df['score'].sort_values().drop_duplicates(keep='first')
        hist, bin = np.histogram(score.values, bins=10)
        return bin[0]

    def cut_edges(self, df, threshold):
        if df[df['score'] > threshold].empty:
            return df

        ridx = []
        N = df.shape[0]
        for i in range(N // 2):
            if df.iloc[i]['score'] <= threshold:
                ridx.append(df.index[i])
            if df.iloc[N - i - 1]['score'] <= threshold:
                ridx.append(df.index[N - i - 1])
        return df.drop(ridx)

    def merge_tags(self, df, distance=20, threshold=None):
        if df.shape[0] < 2:
            start, end = df.iloc[0]['start'], df.iloc[0]['end']
            width = end - start
            df_res = pd.DataFrame(data=[['[{},{})'.format(start, end), width, df.iloc[0]['score']]], columns=['region', 'width', 'score'])
            return df_res
        diff = df['start'].diff()
        idx = diff[diff > distance].index
        N = len(idx)
        if N == 0:
            start, end = df['start'].min(), df['end'].max()
            width = end - start
            df_res = pd.DataFrame(data=[['[{},{})'.format(start, end), width, df['score'].mean()]],
                         columns=['region', 'width', 'score'])
            return df_res
        indeces = np.zeros(N + 2).astype(int)
        indeces[0] = df.index[0]
        indeces[-1] = df.index[-1] + 1
        indeces[1:-1] = idx

        dfs = []
        for i in range(N + 1):
            sidx = indeces[i]
            eidx = indeces[i + 1] - 1
            df_sub = df.loc[sidx: eidx, :]
            if df_sub.shape[0] > 2:
                df_sub = self.cut_edges(df_sub, threshold)
            start = df_sub['start'].min()
            end = df_sub['end'].max()
            score = df_sub['score'].mean()
            width = end - start
            df_new = pd.DataFrame(data=[['[{},{})'.format(start, end), width, score]], columns=['region', 'width', 'score'])
            dfs.append(df_new)
        df_res = pd.concat(dfs)
        return df_res

    def merge_versions(self):
        df = pd.read_excel('correlation_{}_{}_v0.xlsx'.format(self.band, self.csize), index_col=0)

        for i in range(3):
            df_new = pd.read_excel('correlation_{}_{}_v{}.xlsx'.format(self.band, self.csize, i + 1), index_col=0)
            df.loc[df_new.index, 'tss v{}'.format(i + 1)] = df_new['tss']
            df.loc[df_new.index, 'spec v{}'.format(i + 1)] = df_new['spec']
            df.loc[df_new.index, 'corr (spearman) v{}'.format(i + 1)] = df_new['corr (spearman)']
        df.to_excel('correlation_{}_{}.xlsx'.format(self.band, self.csize))

    def run(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con = sqlite3.connect(fpath)
        df_ref = pd.read_sql_query("SELECT chromosome, start, end, strand, gene_name FROM 'gencode.v30lift37' WHERE "
                                   "feature='transcript'", con, index_col='gene_name')

        fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con_rna = sqlite3.connect(fpath_rna)

        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue_{}_{}_v{}.db'.format(self.band, self.csize, self.version))
        con_fan = sqlite3.connect(fpath_fan)
        tlist_fan = Database.load_tableList(con_fan)

        # fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_out.db')
        # con_rna = sqlite3.connect(fpath_rna)

        # fpath_fan = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_{}_v{}.db'.format(self.band, self.csize, self.version))
        # con_fan = sqlite3.connect(fpath_fan)
        # tlist_fan = Database.load_tableList(con_fan)

        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'correlation_{}_{}_v{}.db'.format(self.band, self.csize, self.version))
        out_con = sqlite3.connect(out_path)

        for tissue, tname_rna in tlist_fan.items():
            print(tissue)
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con_fan, index_col='gene_name')
            for rtname in tname_rna:
                df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(rtname), con_rna, index_col='gene_name')

                common_genes = sorted(list(set.intersection(set(df_fan.index), set(df_rna.index))))
                df_res = pd.concat([df_ref.loc[common_genes], df_fan.loc[common_genes, 'score'],
                                    df_fan.loc[common_genes, 'scores'], df_rna.loc[common_genes, 'FPKM'],
                                    df_fan.loc[common_genes, 'widths'], df_fan.loc[common_genes, 'tss']], axis=1)
                df_res.to_sql(rtname, out_con, if_exists='replace')

    def corr_versions(self):
        for band in [100, 500]:
            for cluster_size in [20, 100]:
                fname = 'correlation_{}_{}.xlsx'.format(band, cluster_size)
                df = pd.read_excel(fname)
                df = df.dropna(subset=['corr (spearman)'])

                for i in range(1, 4):
                    df = df.dropna(subset=['corr (spearman) v{}'.format(i)])
                    ccoeff = pearsonr(df['corr (spearman)'], df['corr (spearman) v{}'.format(i)])
                    print('[v0 & v{}] fname: {}, coeff: {:0.2f}'.format(i, fname, ccoeff[0]))

    def processInput_get_vector(self, df_sub, gene):
        return [gene, ','.join(df_sub['score'].astype(str)), ','.join(df_sub['start'].astype(str))]

    def get_vector(self, fname):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', fname)
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        out_path = fpath.replace('.db', '.xlsx')
        writer = pd.ExcelWriter(out_path, engine='xlsxwriter')

        for tname in tlist:
            print(tname)
            if "/" in tname:
                continue
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            df_grp = df.groupby('gene_name')
            print(tname, len(df_grp))

            contents = Parallel(n_jobs=self.num_cores)(delayed(self.processInput_get_vector)(df_sub, gene) for gene, df_sub in df_grp)
            df_res = pd.DataFrame(data=contents, columns=['gene', 'scores', 'starts'])
            df_res.to_excel(writer, sheet_name=tname, index=None)
        writer.save()
        writer.close()

    def processInput_correlation_mir_gene0(self, df_mir, mir, i):
        if i % 100 == 0 or i + 1 == df_mir.shape[0]:
            print('{:0.2f}%'.format((i + 1) * 100 / df_mir.shape[0]))

        score = df_mir.loc[mir, 'scores']
        if isinstance(score, str):
            scores_mir = np.array(list(map(int, score.split(','))))
        else:
            scores_mir = np.array(score)
        return scores_mir.mean()

    def processInput_correlation_mir_gene1(self, df_gene, gene, i):
        if i % 1000 == 0 or i + 1 == df_gene.shape[0]:
            print('{:0.2f}%'.format((i + 1) * 100 / df_gene.shape[0]))

        score = df_gene.loc[gene, 'scores']
        if isinstance(score, str):
            scores_gene = np.array(list(map(int, score.split(','))))
        else:
            scores_gene = np.array(score)
        return scores_gene.mean()

    def processInput_correlation_mir_gene2(self, mmean, df_res_gene, gene):
        gmean = df_res_gene.loc[gene]
        ccoeff = pearsonr(mmean, gmean)
        return [gene, ccoeff[0], ccoeff[1]]
        # return [gene, ccoeff.correlation]

    def get_high_correlated_genes(self, band, scope):
        fpath = 'correlation_{}_{}.xlsx'.format(band, scope)
        df = pd.read_excel(fpath, index_col=0)
        return df[df['corr (spearman)'] > 0.6].to_excel(fpath.replace('.xlsx', '_high.xlsx'))

    def get_scores_by_tissues(self, band):
        fpath_mir = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_mir_{}.xlsx'.format(band))
        df_mirs = pd.ExcelFile(fpath_mir)
        tissues = df_mirs.sheet_names

        # high_correlated_genes = pd.read_excel('correlation_{}_20.xlsx'.format(band), index_col=0).index
        fpath_gene = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}.xlsx'.format(band))
        df_genes = pd.ExcelFile(fpath_gene)

        mirs, genes = [], []
        for tissue in tissues:
            df_mir = df_mirs.parse(tissue, index_col=0)
            df_gene = df_genes.parse(tissue, index_col=0)
            mirs.append(set(df_mir.index))
            genes.append(set(df_gene.index))
        
        mirs = sorted(list(set.intersection(*mirs)))
        genes = sorted(list(set.intersection(*genes)))

        print(len(mirs))
        print(len(genes))

        df_res_mir, df_res_gene = [], []
        for tissue in tissues:
            print(tissue)
            df_mir = df_mirs.parse(tissue, index_col=0).loc[mirs]
            df_gene = df_genes.parse(tissue, index_col=0).loc[genes]

            scores_mir_mean = Parallel(n_jobs=self.num_cores)(delayed(self.processInput_correlation_mir_gene0)(df_mir, mir, i) for i, mir in enumerate(df_mir.index))
            df_res_mir.append(pd.Series(scores_mir_mean, index=df_mir.index))
            print('mir is done')
            scores_gene_mean = Parallel(n_jobs=self.num_cores)(delayed(self.processInput_correlation_mir_gene1)(df_gene, gene, i) for i, gene in enumerate(df_gene.index))
            df_res_gene.append(pd.Series(scores_gene_mean, index=df_gene.index))
            print('gene is done')

        df_res_mir = pd.concat(df_res_mir, axis=1)
        df_res_mir.columns = tissues
        df_res_gene = pd.concat(df_res_gene, axis=1)
        df_res_gene.columns = tissues

        df_res_mir.to_excel(fpath_mir.replace('.xlsx', '_score_vector.xlsx'))
        df_res_gene.to_excel(fpath_gene.replace('.xlsx', '_score_vector.xlsx'))

    def correlation_mir_gene(self, band):
        fpath_mir = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_mir_{}_score_vector.xlsx'.format(band))
        df_res_mir = pd.read_excel(fpath_mir, index_col=0)

        fpath_gene = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector.xlsx'.format(band))
        df_res_gene = pd.read_excel(fpath_gene, index_col=0)

        fpath_out = fpath_gene.replace('_vector.xlsx', '_corr.xlsx')
        writer_corr = pd.ExcelWriter(fpath_out, engine='xlsxwriter')

        fpath_out = fpath_gene.replace('_vector.xlsx', '_pval.xlsx')
        writer_pval = pd.ExcelWriter(fpath_out, engine='xlsxwriter')

        dfs_corr = []
        dfs_pval = []
        N = df_res_mir.shape[0]
        for i, mir in enumerate(df_res_mir.index):
            print('{} [{} / {}]'.format(mir, i + 1, N))
            mmean = df_res_mir.loc[mir]

            table = Parallel(n_jobs=self.num_cores)(delayed(self.processInput_correlation_mir_gene2)(mmean, df_res_gene, gene) for gene in df_res_gene.index)
            df = pd.DataFrame(data=table, columns=['gene', 'ccoeff', 'p-value']).set_index('gene')
            df.columns = [mir, mir]
            print(df.iloc[:, 0])
            print(df.iloc[:, 1])
            dfs_corr.append(df.iloc[:, 0])
            dfs_pval.append(df.iloc[:, 1])

        for df, writer in zip([dfs_corr, dfs_pval], [writer_corr, writer_pval]):
            df = pd.concat(df, axis=1)
            df.to_excel(writer)
            writer.save()
            writer.close()

    def intersection_versions(self):
        df = pd.read_excel('correlation_v0.xlsx', index_col=0)
        df = df[df['corr (spearman)'] > 0.6]
        genes0 = set(df.index)

        genes = []
        num = []
        for i in range(1, 4):
            df = pd.read_excel('correlation_v{}.xlsx'.format(i), index_col=0)
            df = df[df['corr (spearman)'] > 0.6]
            genes.append(set(df.index))
            num.append(df.shape[0])

        N = len(genes0)
        for i, gene in enumerate(genes):
            common_genes = set.intersection(genes0, gene)
            print('(v0) vs. (v{}) {:0.2f}%'.format(i + 1, 100 * len(common_genes) / N))

    def figure(self):
        import matplotlib.pyplot as plt

        genes = []
        dfs = []
        for i in range(3):
            df = pd.read_excel('correlation_v{}.xlsx'.format(i + 1), index_col=0)
            df = df.dropna(subset=['corr (spearman)'])
            dfs.append(df)
            genes.append(set(df.index))

        genes = sorted(list(set.intersection(*genes)))
        corr1 = dfs[0].loc[genes, 'corr (spearman)'].values
        corr2 = dfs[2].loc[genes, 'corr (spearman)'].values
        ccoeff = pearsonr(corr1, corr2)
        print('correlation coeff: {:0.4f}'.format(ccoeff[0]))

        corr_diff = (dfs[0].loc[genes, 'corr (spearman)'] - dfs[2].loc[genes, 'corr (spearman)']) / dfs[2].loc[genes, 'corr (spearman)']
        genes = corr_diff[corr_diff.abs() > 0.5].index

        N = 6
        for i, gene in enumerate(genes):
            fpkm = list(map(float, dfs[0].loc[gene, 'fpkms'].split(',')))
            scores1 = list(map(float, dfs[0].loc[gene, 'scores'].split(':')))
            scores3 = list(map(float, dfs[2].loc[gene, 'scores'].split(':')))

            xaxis = range(len(scores1))

            j = i % N
            if j == 0:
                plt.figure(figsize=(12, 8))
            fname = 'cluster{:02d}.png'.format(int(i // N))
            plt.subplot(2, 3, j + 1)
            plt.plot(xaxis, scores1, label='cluster (v1)')
            plt.plot(xaxis, scores3, label='max (v3)')
            plt.plot(xaxis, fpkm, label='fpkm')
            plt.title(gene)
            plt.legend()
            if j == N - 1:
                plt.savefig(os.path.join(self.root, 'database/Fantom/v5/tissues/out/figures', fname))
                plt.close()


if __name__ == '__main__':
    cor = Correlation2()
    if cor.hostname == 'mingyu-Precision-Tower-7810':
        cor.to_server()

    else:
        for band in [500]:
            cor.band = band
            # cor.get_vector('human_cell_line_hCAGE_{}.db'.format(band))
            # cor.get_vector('human_cell_line_hCAGE_mir_{}.db'.format(band))

            for cluster_size in [100]:
                cor.csize = cluster_size
                for i in range(0):
                    cor.version = i
                    cor.fantom_unique_gene()
                    cor.run()
                    cor.correlation()
                cor.merge_versions()
                cor.get_high_correlated_genes(band, cluster_size)

            cor.get_scores_by_tissues(band)
            cor.correlation_mir_gene(band)
