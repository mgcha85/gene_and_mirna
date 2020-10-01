import sqlite3
import os
import socket
import pandas as pd
import numpy as np


class comparison:
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

    def wilcox(self):
        from scipy.stats import wilcoxon
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()

        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out')
        con_out = sqlite3.connect(os.path.join(dirname, 'regression_100_nz_wilcox.db'))

        con = sqlite3.connect(os.path.join(dirname, 'regression_100_nz.db'))
        df = pd.read_sql("SELECT * FROM 'corr'", con, index_col='transcript_id')

        # con_inter = sqlite3.connect(os.path.join(dirname, 'regression_100_nz_inter.db'))
        # df_inter = pd.read_sql("SELECT * FROM 'corr'", con_inter, index_col='transcript_id')

        # genes = set.intersection(set(df.index), set(df_inter.index))
        # df_inter = df_inter.loc[genes]
        # df = df.loc[genes]
        #
        # diff = np.abs(np.subtract(df, df_inter))
        # df_diff = pd.DataFrame(data=diff, columns=df.index, index=df.columns)

        report = []
        for i in range(df.shape[1]):
            mir1 = df.columns[i]
            for j in range(df.shape[1]):
                if i == j:
                    continue

                mir2 = df.columns[j]
                x = df.iloc[:, i]
                y = df.iloc[:, j]
                diff = (x - y).sum()

                if diff != 0:
                    T, p_value = wilcoxon(x, y, alternative='greater')
                    if diff > 0 and p_value < 0.01:
                        report.append([mir1, mir2, T, p_value, diff])
                    else:
                        report.append([mir1, mir2, None, None, 0])
                else:
                    report.append([mir1, mir2, None, None, 0])

            df_rep = pd.DataFrame(report, columns=['miRNA1', 'miRNA2', 'sum', 'p-value', 'diff'])
            df_rep.sort_values(by=['diff']).to_sql('wilcox_corr_lasso', con_out, index=None, if_exists='replace')
        # Parallel(n_jobs=num_cores)(delayed(processInput)(fname) for fname in flist)

    def convert_pre_to_mir(self):
        fpath_ref = os.path.join(self.root, 'database/mirbase/22', 'mirbase.db')
        con_ref = sqlite3.connect(fpath_ref)
        df_ref = pd.read_sql("SELECT * FROM 'hsa_liftover_hg19'", con_ref)
        df_ref = df_ref.drop_duplicates(subset=['Name'])

        df_ref['Name'] = df_ref['Name'].str.lower()
        df_ref = df_ref.set_index('Name', drop=True)

        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)

        for tname in ['miRTartBase_hsa', 'target_scan_grp']:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='miRNA')
            for mir in df.index:
                if mir != 'hsa-mir-1244':
                    continue
                if mir not in df_ref.index:
                    continue
                df_mir = df_ref.loc[mir]
                if isinstance(df_mir, pd.DataFrame):
                    df_mir = df_mir.reset_index()
                    idx = df_mir[df_mir['feature'] == 'miRNA_primary_transcript'].index[0]
                    premir = df_mir.loc[idx, 'Name']
                else:
                    if df_mir['feature'] == 'miRNA_primary_transcript' or len(mir.split('-')) < 4:
                        premir = mir
                    else:
                        premir = '-'.join(mir.split('-')[:-1])
                df.loc[mir, 'pre-miRNA'] = premir
            df.to_sql(tname, con, if_exists='replace')

    def get_ref_transcripts(self):
        from Database import Database

        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)

        fpath_ref = os.path.join(self.root, 'database/gencode', 'gencode.v32lift37.annotation_attr_spt.db')
        con_ref = sqlite3.connect(fpath_ref)
        dfs = []
        for tname in Database.load_tableList(con_ref):
            chr, str = tname.split('_')
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con_ref)
            df['chromosome'] = chr
            df['strand'] = str
            dfs.append(df)

        con_mir = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz.db'))
        mir_list = pd.read_sql("SELECT miRNA FROM 'Y'", con_mir)
        df_ref = pd.concat(dfs)
        for tname in ['miRTartBase_hsa', 'target_scan_grp']:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='pre-miRNA')
            mir = set.intersection(set(mir_list['miRNA']), set(df.index))
            df = df.loc[mir]
            df = df.reset_index()
            contents = []

            for idx, mir in enumerate(df.index):
                if (idx + 1) % 100 == 0 or idx + 1 == df.shape[0]:
                    print('{} / {}'.format(idx + 1, df.shape[0]))
                genes = df.loc[mir, 'genes'].split(';')
                for gene in genes:
                    contents.append(df_ref[df_ref['gene_name'] == gene])
            df_res = pd.concat(contents)

            con_out = sqlite3.connect(os.path.join(self.root, 'database', '{}_tr.db'.format('_'.join(tname.split('_')[:-1]))))
            for str, df_str in df_res.groupby('strand'):
                for chr, df_chr in df_str.groupby('chromosome'):
                    df_chr.drop(['strand', 'chromosome'], axis=1).to_sql('{}_{}'.format(chr, str), con_out, index=False, if_exists='replace')

    def compare_corrlation_pair(self):
        from scipy.stats import mannwhitneyu
        from itertools import combinations
        from Database import Database

        fpath = os.path.join(self.root, "database/Fantom/v5/cell_lines/out", 'regression_100_nz.db')
        con = sqlite3.connect(fpath)
        pairs = [('corr', 'corr_ts'), ('corr', 'corr_mt')]

        fpath = os.path.join(self.root, "database", 'target_genes.db')
        con_other = sqlite3.connect(fpath)

        fpath_mt = os.path.join(self.root, "database", 'miRTartBase_tr.db')
        con_mt = sqlite3.connect(fpath_mt)
        dfs = []
        for tname in Database.load_tableList(con_mt):
            dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con_mt))
        df_mt_tr = pd.concat(dfs).drop_duplicates(subset=['transcript_id'])

        df_ts = pd.read_sql("SELECT * FROM 'target_scan_grp'", con_other).dropna(subset=['pre-miRNA'])
        df_mt = pd.read_sql("SELECT * FROM 'miRTartBase_hsa'", con_other).dropna(subset=['pre-miRNA'])
        df_la = pd.read_sql("SELECT * FROM 'result'", con, index_col='miRNA')

        dfs_corr = {'corr': pd.read_sql("SELECT * FROM 'corr'", con, index_col='transcript_id'),
                'corr_ts': pd.read_sql("SELECT * FROM 'corr_ts'", con, index_col='transcript_id'),
                'corr_mt': pd.read_sql("SELECT * FROM 'corr_mt'", con, index_col='transcript_id')}
        for key, df_corr in dfs_corr.items():
            dfs_corr[key] = df_corr[~df_corr.index.duplicated(keep='first')]

        dfs_res = []
        for df_other, label in zip([df_ts, df_mt], ['ts', 'mi']):
            cmirs = set.intersection(set(df_other['pre-miRNA']), set(df_la.index))

            result = []
            data_other = {}
            data_la = {}
            for i, mir in enumerate(cmirs):
                print('{} / {}'.format(i+1, len(cmirs)))

                data_other[mir] = []
                idx = df_other[df_other['pre-miRNA'] == mir].index
                for i in idx:
                    genes = df_other.loc[i, 'genes'].split(';')
                    for gene in genes:
                        for tr in df_mt_tr[df_mt_tr['gene_name'] == gene]['transcript_id']:
                            data_other[mir].append(dfs_corr['corr_mt'].loc[tr, mir])
                # Lasso
                data_la[mir] = []
                for tr in df_la.loc[mir, 'Transcripts'].split(';'):
                    data_la[mir].append(dfs_corr['corr'].loc[tr, mir])

                stat, p = mannwhitneyu(np.array(data_other[mir]), np.array(data_la[mir]))
                result.append([mir, stat, p, len(data_other[mir]), len(data_la[mir])])
            dfs_res.append(pd.DataFrame(result, columns=['miRNA', 'stat ({})'.format(label), 'p-value ({})'.format(label), '# genes ({})'.format(label), '# genes (Lasso)']).set_index('miRNA'))
        df_res = pd.concat(dfs_res, axis=1)
        df_res.index.name = 'miRNA'
        df_res.to_excel(os.path.join(self.root, "database/Fantom/v5/cell_lines/out", 'mannwhitneyu.xlsx'))

    def run(self):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz_union.db')
        con_union = sqlite3.connect(fpath)

        df_cor_union = pd.read_sql("SELECT * FROM 'corr'", con_union, index_col='transcript_id')

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz_inter.db')
        con_inter = sqlite3.connect(fpath)

        df_cor_inter = pd.read_sql("SELECT * FROM 'corr'", con_inter, index_col='transcript_id')

        plt.subplot(211)
        xaxis = range(df_cor_union.shape[1])
        plt.plot(xaxis, df_cor_union.mean())
        plt.xticks(xaxis[::10], df_cor_union.columns[::10], fontsize=6, rotation=30)
        plt.title('corr coeff by a miRNA [union]')
        plt.grid()

        plt.subplot(212)
        xaxis = range(df_cor_inter.shape[1])
        plt.plot(xaxis, df_cor_inter.mean())
        plt.xticks(xaxis[::10], df_cor_inter.columns[::10], fontsize=6, rotation=30)
        plt.title('corr coeff by a miRNA [intersection]')
        plt.grid()
        plt.show()

        df_coef_union = pd.read_sql("SELECT * FROM 'coefficient'", con_union, index_col='transcript_id')
        df_coef_inter = pd.read_sql("SELECT * FROM 'coefficient'", con_inter, index_col='transcript_id')
        z_count_union = (df_coef_union != 0).astype(int).sum(axis=0)
        z_count_inter = (df_coef_inter != 0).astype(int).sum(axis=0)
        zratio_union = 100 * z_count_union / df_coef_union.shape[0]
        zratio_inter = 100 * z_count_inter / df_coef_inter.shape[0]
        print('zratio_inter: {:0.2f}%, zratio_union: {:0.2f}%'.format(zratio_inter.mean(), zratio_union.mean()))
        print('inter mean: {:0.2f}, union mean: {:0.2f}'.format(df_cor_inter.mean().mean(), df_cor_union.mean().mean()))


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

    comp = comparison(root)
    if hostname == 'DESKTOP-DLOOJR6' or hostname == '-1NLOLK4':
        comp.to_server(root, "")
    else:
        # comp.convert_pre_to_mir()
        # comp.get_ref_transcripts()
        comp.compare_corrlation_pair()