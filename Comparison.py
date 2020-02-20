import sqlite3
import os
import socket
import pandas as pd
import numpy as np


class comparison:
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

        which = 'stokes'
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
    comp = comparison()
    if comp.hostname == 'mingyu-Precision-Tower-7810':
        comp.to_server()
    else:
        comp.wilcox()
