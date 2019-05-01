import os
import sqlite3
import pandas as pd
import socket


class Fantom_RNA:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.cells = []

    def tpm(self, df):
        print('do something')

    def run(self):
        fpath_out = os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc_full_out.db')
        con_out = sqlite3.connect(fpath_out)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc_full.db')
        con = sqlite3.connect(fpath)

        fpath = os.path.join(self.root, 'database/RNA-seq', 'rna_seq.db')
        con_rna = sqlite3.connect(fpath)
        for cell in self.cells:
            df = pd.read_sql_query("SELECT * FROM '{}'".format(cell), con)
            for idx in df.index:
                chromosome = df.loc[idx, 'chromosome']
                start = df.loc[idx, 'start']
                end = df.loc[idx, 'end']
                strand = df.loc[idx, 'strand']

                df_rna = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND strand='{}' AND start>{end} "
                                           "AND end<{start}".format(cell, chromosome, strand, start=start, end=end),
                                           con_rna)
                df.loc[idx, 'TPM'] = self.tpm(df_rna)

            df.to_sql(cell, con_out, index=None, if_exists='replace')

