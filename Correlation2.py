import pandas as pd
import sqlite3
import os
import socket
from Database import Database
import numpy as np


class Correlation2:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.cells = []

    def pre_data(self):
        fpath = os.path.join(self.root, 'software/stringtie-1.3.6', 'genes.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        out_con = sqlite3.connect(fpath.replace('.db', '_2.db'))

        for tname in tlist:
            df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            attribute = df_rna['attribute'].str.split('; ')

            pkm = []
            for attr in attribute:
                dict = {}
                for a in attr[:-1]:
                    key, value = a.split(' ')
                    dict[key] = value.replace('"', '')

                if 'gene_name' not in dict:
                    pkm.append(np.nan)
                    continue
                pkm.append(dict['gene_name'])

            df_res = pd.DataFrame(data=pkm, columns=['gene_name'])
            df_rna = pd.concat([df_rna, df_res], axis=1)
            df_rna.drop(['attribute', 'frame'], axis=1).to_sql(tname, out_con, index=None, if_exists='append')

    def run(self):
        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        fpath_ref = os.path.join(self.root, 'software/stringtie-1.3.6', 'genes_2.db')
        con_ref = sqlite3.connect(fpath_ref)

        for tname in tlist:
            df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)

            for idx in df_rna.index:
                chromosome = df_rna.loc[idx, 'chromosome']
                strand = df_rna.loc[idx, 'strand']
                gene_name = df_rna.loc[idx, 'gene_name']

                gene_tname = 'genes_{}_{}'.format(chromosome, strand)
                df_ref = pd.read_sql_query("SELECT * FROM '{}' WHERE source='protein_coding' AND feature='transcript' "
                                           "AND gene_name={}".format(gene_tname, gene_name), con_ref)


if __name__ == '__main__':
    cor = Correlation2()
    cor.run()
