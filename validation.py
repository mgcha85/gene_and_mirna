import sqlite3
import os
import pandas as pd
import socket
from scipy.stats import wilcoxon
from Database import Database
import numpy as np


class validation:
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

    def comparison_by_size(self):
        dfList = {}
        for pair in [[100, 300], [300, 500], [100, 500]]:
            dfs = []
            for i, s in enumerate(pair):
                fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}_nz.db'.format(s))
                con = sqlite3.connect(fpath)
                df = pd.read_sql("SELECT miRNA, gene_name AS 'gene_name_{}' FROM 'result'".format(s), con, index_col='miRNA')
                df['gene_name_{}'.format(s)] = df['gene_name_{}'.format(s)].str.split(';')
                dfs.append(df)

            df = pd.concat(dfs, axis=1)
            df = df.dropna(how='any')
            df_res = pd.DataFrame(index=df.index)

            stats = pd.Series(index=range(df.shape[0]))
            for i, mir in enumerate(df.index):
                sizes = []
                dfs = []
                for col in df.columns:
                    genes = set(df.loc[mir, col])
                    sizes.append(len(genes))
                    dfs.append(genes)

                intersection = set.intersection(*dfs)
                union = set.union(*dfs)

                min_idx = sizes.index(min(sizes))
                max_idx = sizes.index(max(sizes))
                stats[i] = sizes[max_idx]
                df_res.loc[mir, '# {}'.format(pair[0])] = sizes[0]
                df_res.loc[mir, '# {}'.format(pair[1])] = sizes[1]
                df_res.loc[mir, '{} ∩ {}'.format(pair[0], pair[1])] = len(intersection)
                df_res.loc[mir, '∩/# smaller'] = len(intersection) / sizes[min_idx]
                df_res.loc[mir, '∩/# larger'] = len(intersection) / sizes[max_idx]
                df_res['larger percentage'] = df_res[['∩/# smaller', '∩/# larger']].max(axis=1)
            dfList['{} vs {}'.format(*pair)] = df_res
            print('[{} vs {}] mean: {:0.4f}, median: {:0.4f}, min: {:0.4f}, max: {:0.4f}, q20: {:0.4f}, q40: {:0.4f}, q60: {:0.4f}'.format(*pair, df_res['larger percentage'].mean(), df_res['larger percentage'].median(), df_res['larger percentage'].min(), df_res['larger percentage'].max(), df_res['larger percentage'].quantile(q=0.2), df_res['larger percentage'].quantile(q=0.4), df_res['larger percentage'].quantile(q=0.6)))
        writer = pd.ExcelWriter(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'comparison_by_size.xlsx'), engine='xlsxwriter')
        for key, df in dfList.items():
            df.to_excel(writer, sheet_name=key)
            df.to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'comparison_by_size.xlsx'))
        writer.save()
        writer.close()

    def wilcox(self, opt):
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()

        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', opt)
        flist = os.listdir(dirname)
        flist = sorted([x for x in flist if x.endswith('.db') and 'test' in x])

        def processInput(fname):
            con_out = sqlite3.connect(os.path.join(dirname, 'regression_100_cv.db'))
            fpath = os.path.join(dirname, fname)
            con = sqlite3.connect(fpath)
            df = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='tid')
            df_x = pd.read_sql("SELECT * FROM 'X'", con, index_col='tid').T
            df_y = pd.read_sql("SELECT * FROM 'Y'", con, index_col='miRNA').T

            y = np.matmul(df_x.values, df.values)
            diff = np.abs(np.subtract(df_y.values, y))
            df_diff = pd.DataFrame(data=diff, columns=df_y.columns, index=df_y.index)

            report = []
            for i in range(df.shape[1]):
                mir1 = df_diff.columns[i]
                for j in range(df.shape[1]):
                    if i == j:
                        continue

                    mir2 = df_diff.columns[j]
                    x = df_diff.iloc[:, i]
                    y = df_diff.iloc[:, j]
                    diff = (x - y).sum()

                    if diff != 0:
                        T, p_value = wilcoxon(x, y, alternative='greater')
                        if diff > 0 and p_value < 0.01:
                            report.append([mir1, mir2, T, p_value, diff])
                        elif diff < 0:
                            T, p_value = wilcoxon(x, y, alternative='less')
                            if p_value < 0.01:
                                report.append([mir1, mir2, T, p_value, diff])
                            else:
                                report.append([mir1, mir2, None, None, 0])
                        else:
                            report.append([mir1, mir2, None, None, 0])
                    else:
                        report.append([mir1, mir2, None, None, 0])
            df_rep = pd.DataFrame(report, columns=['miRNA1', 'miRNA2', 'sum', 'p-value', 'diff'])
            df_rep.sort_values(by=['diff']).to_sql(os.path.splitext(fname)[0], con_out, index=None, if_exists='replace')

        for fname in flist:
            processInput(fname)
        # Parallel(n_jobs=num_cores)(delayed(processInput)(fname) for fname in flist)

    def find_identical_distribution(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', 'regression_100_cv.db')
        con = sqlite3.connect(fpath)

        contents = []
        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_cv_id.xlsx')
        writer = pd.ExcelWriter(out_path, engine='xlsxwriter')

        for tname in Database.load_tableList(con):
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='index')
            ridx, cidx = np.where(df.values < 0.05)
            print('{:0.2f}%'.format(100 * len(ridx) / df.size))
            for r, c in zip(ridx, cidx):
                mir1 = df.columns[c]
                mir2 = df.index[r]
                p_val = df.iloc[r, c]
                contents.append([mir1, mir2, p_val])
            df_res = pd.DataFrame(contents, columns=['miRNA1', 'miRNA2', 'p-value'])
            df_res.to_excel(writer, sheet_name=tname, index=None)
        writer.save()
        writer.close()

    def is_similar(self, opt):
        con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', opt, 'regression_100_cv.db'))
        tlist = sorted(Database.load_tableList(con))

        dfs = []
        lengths = np.zeros(len(tlist))
        for i, tname in enumerate(tlist):
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            df = df[df['p-value'] < 0.05]
            df['pair'] = df['miRNA1'] + ';' + df['miRNA2']
            df = df.set_index('pair')
            dfs.append(df)
            lengths[i] = df.shape[0]

        n = len(dfs)
        report = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue

                pdf = dfs[i]
                df = dfs[j]

                inter = set.intersection(set(pdf.index), set(df.index))
                report[i, j] = (100 * len(inter) / lengths.mean())
        df_res = pd.DataFrame(report)
        df_res.to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', opt, 'regression_100_cv.xlsx'), index=None)

    def quick_sort(self, df):
        if df.shape[0] <= 1:
            return df

        n = int(df.shape[0] * .75)
        frow = df.iloc[n:n+1, :]
        pivot = frow.iloc[0]['miRNA1']

        df_sub = df.drop(df.index[n])
        df_sub = df_sub[df_sub['miRNA1'] == pivot]

        g0 = df_sub[df_sub['diff'] == 0].index
        g1 = df_sub[df_sub['diff'] > 0].index
        g2 = df_sub[df_sub['diff'] < 0].index
        return pd.concat([self.quick_sort(df.loc[g1]), frow, df.loc[g0], self.quick_sort(df.loc[g2])])

    def sort_diff(self, opt):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', opt)
        con = sqlite3.connect(os.path.join(dirname, 'regression_100_cv.db'))
        con_out = sqlite3.connect(os.path.join(dirname, 'regression_100_cv_sorted.db'))

        summary = []
        for tname in Database.load_tableList(con):
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            # df = pd.read_sql("SELECT * FROM '{}' WHERE miRNA1='{}' ORDER BY diff DESC".format(tname, mir), con)
            df = self.quick_sort(df)
            df.to_sql(tname, con_out, if_exists='replace', index=None)
            # summary.append([*df[df['diff'] > 0]['miRNA2'], df['miRNA1'].iloc[0],  *df[df['diff'] == 0]['miRNA2'], *df[df['diff'] < 0]['miRNA2']])
            pivot = df.iloc[0]['miRNA1']
            summary.append([*df[df['diff'] > 0]['miRNA2'], pivot, *df[df['diff'] == 0]['miRNA2'], *df[df['diff'] < 0]['miRNA2']])
        df_sum = pd.DataFrame(data=summary)
        df_sum.T.to_excel('summary.xlsx', index=None)


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

    val = validation(root)
    if hostname == 'DESKTOP-DLOOJR6' or hostname == '-1NLOLK4':
        val.to_server(root, "")
    else:
        val.comparison_by_size()
        # for opt in ['nz']:
            # val.wilcox(opt)
            # val.is_similar(opt)
            # val.sort_diff(opt)
        # val.find_identical_distribution()
