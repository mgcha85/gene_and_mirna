import pandas as pd
import sqlite3
import os
import socket
from Database import Database
import numpy as np
from joblib import Parallel, delayed
import multiprocessing


class Tss_map:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        # self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293']
        self.cell_lines = self.get_cell_lines()[300:]
        self.chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                      'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
        # 26

    def get_cell_lines(self):
        from Database import Database
        return Database.load_tableList(self.connect_cage_tags())

    def connect_cage_tags(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc_full.db')
        return sqlite3.connect(fpath, check_same_thread=False)

    def load_all_tags(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc_loc.csv')
        return pd.read_csv(fpath)

    def run(self):
        df__ = self.load_all_tags()
        df_grp__ = df__.groupby('chromosome')

        fpath = os.path.join(self.root, 'Papers/Tss_map', 'Tss_map_table.db')
        con = sqlite3.connect(fpath)

        out_path = os.path.join(self.root, 'Papers/Tss_map', 'Tss_map.db')
        con_out = sqlite3.connect(out_path)
        for chrom in self.chrom:
            df_sub__ = df_grp__.get_group(chrom)

            for strand in ['+', '-']:
                print(chrom, strand)
                df_sub = df_sub__[df_sub__['strand'] == strand]
                df_sub = df_sub.sort_values('start')
                df_sub.index = df_sub['start'].astype(str) + ';' + df_sub['end'].astype(str)

                for cline in self.cell_lines:
                    df = pd.read_sql("SELECT * FROM '{}_{}_{}'".format(cline, chrom, '+'), con)
                    df.index = df['start'].astype(str) + ';' + df['end'].astype(str)
                    df = df[~df.index.duplicated(keep='first')]
                    df_sub[cline] = df['type']

                tname = '{}_{}'.format(chrom, strand)
                df_sub.to_sql(tname, con_out, if_exists='replace', index=None)


if __name__ == '__main__':
    tm = Tss_map()
    tm.run()
