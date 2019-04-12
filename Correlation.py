import socket
import os
import sqlite3
import pandas as pd
import scipy.stats
import numpy as np


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
        str_map = {'+': 'plus', '-': 'minus'}
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
                                       "strand='{}' AND NOT start>{end} AND NOT end<{start}"
                                       "".format(chromosome, strand, start=start, end=end), con)

            df_gene = pd.read_sql_query("SELECT * FROM 'human_gene' WHERE chromosome='{}' AND NOT "
                                       "start>{end} AND NOT end<{start}"
                                       "".format(chromosome, strand, start=start, end=end), con)
            if not df_mir.empty:
                mirna = df_mir['premiRNA'].str.join(';')
                contents.append(['miRNA', mirna])
            elif not df_gene.empty:
                gene = df_mir['gene'].str.join(';')
                contents.append(['gene', gene])
            else:
                contents.append(['other', None])

        df_type = pd.DataFrame(contents, columns=['tss-type', 'name'])
        df_fantom = pd.concat([df_fantom, df_type], axis=1)
        df_fantom.to_csv(fpath.replace('.csv', '2.csv'), index=None)

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df_fantom = pd.read_csv(fpath)

        df_fantom_mir = df_fantom[df_fantom['tss-type'] == 'miRNA']
        df_fantom_gene = df_fantom[df_fantom['tss-type'] == 'gene']

        matrix = np.zeros((df_fantom_mir.shape[0], df_fantom_gene.shape[0]))
        for i, midx in enumerate(df_fantom_mir.index):
            mrow = df_fantom_mir.loc[midx].iloc[4:]
            for j, gidx in enumerate(df_fantom_gene.index):
                grow = df_fantom_gene.loc[midx].iloc[4:]
                matrix[i, j] = scipy.stats.spearmanr(mrow, grow)
        df_rep = pd.DataFrame(matrix)
        df_rep.to_excel('correlation_report.xlsx')


if __name__ == '__main__':
    cor = Correlation()
    cor.run()