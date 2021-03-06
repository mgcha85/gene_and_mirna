import sqlite3
import os
import pandas as pd
import numpy as np
import socket
import matplotlib.pyplot as plt
import inspect


class Statistic:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def count_lasso_result(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter.xlsx')
        df = pd.read_excel(fpath, index_col=0)
        buffer = np.zeros(df.shape[0])

        for i, mir in enumerate(df.index):
            genes = df.loc[mir, 'GENEs']
            genes = genes.split(';')
            buffer[i] = len(genes)

        bmean = buffer.mean()
        plt.figure(figsize=(10, 8))

        step = 10
        xaxis = np.arange(df.shape[0])
        plt.bar(xaxis, buffer)
        plt.plot(xaxis[[0, -1]], [bmean] * 2, color='r', linestyle=':')
        plt.xticks(xaxis[::step], df.index[::step], rotation=30, fontsize=6)
        plt.grid()
        plt.title('mean: {:0.2f}, std: {:0.2f}, total: {:,.0f}'.format(bmean, buffer.std(), buffer.sum()))

        fname = inspect.currentframe().f_code.co_name
        plt.savefig('/home/mingyu/Pictures/{}.png'.format(fname))
        plt.close()

    def count_gsea_result(self):
        dirname = '/home/mingyu/gsea_home/output/aug28/my_analysis.Gsea.1567011190321'
        flist = os.listdir(dirname)
        flist = [x for x in flist if x.endswith('.xls') and 'CHR' in x]

        n_set = len(flist)
        buffer = np.zeros(n_set)
        gene_set = []
        for i, fname in enumerate(flist):
            df_gsea = pd.read_csv(os.path.join(dirname, fname), sep='\t')
            gset = set(df_gsea['PROBE'].values)
            buffer[i] = len(gset)
            gene_set.append(gset)

        bmean = buffer.mean()
        plt.figure(figsize=(10, 8))

        step = 1
        xaxis = np.arange(n_set)
        plt.bar(xaxis, buffer)
        plt.plot(xaxis[[0, -1]], [bmean] * 2, color='r', linestyle=':')
        plt.xticks(xaxis[::step], flist[::step], rotation=30, fontsize=6)
        plt.grid()
        plt.title('mean: {:0.2f}, std: {:0.2f}, total: {:,.0f}'.format(bmean, buffer.std(), buffer.sum()))

        fname = inspect.currentframe().f_code.co_name
        plt.savefig('/home/mingyu/Pictures/{}.png'.format(fname))
        plt.close()

    def tissue_specific_stats(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out')
        bands = [100, 300, 500]
        df_rep = pd.DataFrame(index=['# transcripts', 'mean', 'std', 'max', 'min', 'median'], columns=['{}bp'.format(x) for x in bands])

        for band in bands:
            fpath = os.path.join(dirname, 'regression_{}.db'.format(band))
            con = sqlite3.connect(fpath)

            writer = pd.ExcelWriter(fpath.replace('.db', '.xlsx'), engine='xlsxwriter')
            for tname in ['coefficient']:
                df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
                sum = (df == 0).sum(axis=0)
                df_rep.loc['# transcripts', '{}bp'.format(band)] = df.shape[0]
                df_rep.loc['mean', '{}bp'.format(band)] = sum.mean()
                df_rep.loc['std', '{}bp'.format(band)] = sum.std()
                df_rep.loc['max', '{}bp'.format(band)] = sum.max()
                df_rep.loc['min', '{}bp'.format(band)] = sum.min()
                df_rep.loc['median', '{}bp'.format(band)] = sum.median()
                df.append(sum, ignore_index=True).to_excel(writer, sheet_name=tname, index=None)
            writer.save()
            writer.close()
        df_rep.to_excel(os.path.join(dirname, 'regression_stats.xlsx'))


if __name__ == '__main__':
    stat = Statistic()
    stat.tissue_specific_stats()
    # stat.count_gsea_result()
    # stat.count_lasso_result()
