import socket
import os
import sqlite3
import pandas as pd
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import pickle as pkl


class Correlation:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def add_type(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df_fantom = pd.read_csv(fpath)

        fpath_fan = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath_fan)

        contents = []
        for idx in df_fantom.index:
            if idx % 1000 == 0 or idx == df_fantom.shape[0] - 1:
                print('{:0.2f}%'.format(100 * (idx + 1) / df_fantom.shape[0]))

            chromosome = df_fantom.loc[idx, 'chromosome']
            strand = df_fantom.loc[idx, 'strand']
            start = df_fantom.loc[idx, 'start']
            end = df_fantom.loc[idx, 'end']

            df_mir = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE chromosome='{}' AND "
                                       "strand='{}' AND tss>={start} AND tss<{end}"
                                       "".format(chromosome, strand, start=start, end=end), con)

            df_gene = pd.read_sql_query("SELECT * FROM 'human_gene' WHERE chromosome='{}' AND start>={start} AND "
                                        "start<{end}".format(chromosome, strand, start=start, end=end), con)
            if not df_mir.empty:
                mirna = ';'.join(df_mir['premiRNA'])
                contents.append(['miRNA', mirna])
            elif not df_gene.empty:
                gene = ';'.join(df_gene['gene'])
                contents.append(['gene', gene])
            else:
                contents.append(['other', None])

        df_type = pd.DataFrame(contents, columns=['tss-type', 'name'])
        df_fantom = pd.concat([df_fantom, df_type], axis=1)
        df_fantom.to_csv(fpath.replace('.csv', '2.csv'), index=None)

    def run(self):
        # fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc2.csv')
        # df_fantom = pd.read_csv(fpath)

        # df_fantom_mir = df_fantom[df_fantom['tss-type'] == 'miRNA']
        # df_fantom_gene = df_fantom[df_fantom['tss-type'] == 'gene']
        # df_fantom_mir.to_excel('miRNA.xlsx', index=None)
        # df_fantom_gene.to_excel('gene.xlsx', index=None)

        df_fantom_mir = pd.read_excel('miRNA.xlsx')
        df_fantom_gene = pd.read_excel('gene.xlsx')

        matrix = np.zeros((df_fantom_mir.shape[0], df_fantom_gene.shape[0]))
        index = []
        columns = []

        for i, gidx in enumerate(df_fantom_gene.index):
            print(i)
            index.append(df_fantom_gene.loc[gidx, 'name'])
            mrow = df_fantom_gene.loc[gidx].iloc[4:-3]
            for j, midx in enumerate(df_fantom_mir.index):
                columns.append(df_fantom_mir.loc[midx, 'name'])
                grow = df_fantom_mir.loc[midx].iloc[4:-3]

                matrix[i, j] = scipy.stats.spearmanr(mrow, grow)[0]
        df_rep = pd.DataFrame(matrix, index=index, columns=columns)
        df_rep.to_excel('correlation_report.xlsx')

    def split(self):
        fname = 'correlation_report'
        df = pd.read_excel('{}.xlsx'.format(fname), index_col=0)
        df.index.name = 'miRNA'

        N = 500  # number of split tables

        out_path = os.path.join(self.root, 'database', '{}.db'.format(fname))
        con = sqlite3.connect(out_path)

        M = int((df.shape[1] + N - 1) / N)
        k = len(str(M))
        for i in range(M):
            idx = i * N
            df_spt = df.iloc[:, idx: idx + N]
            df_spt.to_sql(fname + '_{}'.format(str(i).zfill(k)), con, if_exists='replace')

    def figure(self):
        df = pd.read_excel('correlation_report.xlsx', index_col=0)

        fig = plt.figure()
        # ax = Axes3D(fig)
        ax1 = fig.add_subplot(111, projection='3d')

        x = range(df.shape[1])
        y = range(df.shape[0])
        X, Y = np.meshgrid(x, y)
        values = df.values.flatten()
        Z = np.zeros_like(values)

        # color_values = plt.cm.jet(df.values.tolist())
        # ax1.bar3d(X, Y, df.values, dx=1, dy=1, dz=1, color=color_values)
        ax1.bar3d(X.flatten(), Y.flatten(), Z, dx=1, dy=1, dz=values)
        plt.savefig('correlation_report.png')

    def filter(self):
        from Database import Database

        fpath = os.path.join(self.root, 'database', 'correlation_report.db')
        con = sqlite3.connect(fpath)
        tableList = Database.load_tableList(con)
        thres = 0.8

        dfs = {'GENE': [], 'miRNA': [], 'value': []}
        for tname in tableList:
            print(tname)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con, index_col='miRNA')
            df = df[~df.index.duplicated(keep='first')]

            for idx in df.index:
                for col in df.columns:
                    val = df.loc[idx, col]
                    # if not isinstance(val, float):
                    #     print('see')
                    if abs(val) > thres:
                        dfs['miRNA'].append(idx)
                        dfs['GENE'].append(col)
                        dfs['value'].append(val)
        df_res = pd.DataFrame(dfs)
        df_res.to_excel('correlation_report2_{:.1f}.xlsx'.format(thres), index=None)


if __name__ == '__main__':
    cor = Correlation()
    # cor.add_type()
    # cor.run()
    # cor.split()
    cor.filter()
