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

        fpath = os.path.join(self.root, 'database', 'genecode.db')
        con = sqlite3.connect(fpath)

        for idx in df_fantom.index:
            chromosome = df_fantom.loc[idx, 'chromosome']
            start = df_fantom.loc[idx, 'start']
            end = df_fantom.loc[idx, 'end']

            tname = 'gencode_v28_transcripts_{}_{}'.format(chromosome, strand[df_fantom.loc[idx, 'strand']])
            df_mir = pd.read_sql_query("SELECT * FROM '{}' WHERE 'MIR' IN gene".format(tname), con)
            print(df_mir)
            exit(1)


if __name__ == '__main__':
    cor = Correlation()
    cor.run()
