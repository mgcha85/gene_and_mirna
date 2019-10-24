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


if __name__ == '__main__':
    stat = Statistic()
    stat.count_gsea_result()
    stat.count_lasso_result()
