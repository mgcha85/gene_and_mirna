import socket
import os
import sqlite3
import pandas as pd


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

    def load_fantom(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        return pd.read_csv(fpath)

    def run(self):
        strand = {'+': 'plus', '-': 'minus'}
        df_fantom = self.load_fantom()

        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)

        contents = []
        for idx in df_fantom.index:
            chromosome = df_fantom.loc[idx, 'chromosome']
            start = df_fantom.loc[idx, 'start']
            end = df_fantom.loc[idx, 'end']

            df_mir = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE chromosome='{}' AND "
                                       "strand='{}' AND NOT start>{end} AND NOT end<{start}"
                                       "".format(chromosome, strand, start=start, end=end), con)

            df_gene = pd.read_sql_query("SELECT * FROM 'human_gene' WHERE chromosome='{}' AND strand='{}' AND NOT "
                                       "start>{end} AND NOT end<{start}"
                                       "".format(chromosome, strand, start=start, end=end), con)
            if not df_mir.empty:
                contents.append('miRNA')
            elif not df_gene.empty:
                contents.append('gene')
            else:
                contents.append('other')


if __name__ == '__main__':
    cor = Correlation()
    cor.run()
