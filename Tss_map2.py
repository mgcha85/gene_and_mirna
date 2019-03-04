import pandas as pd
import sqlite3
import os
import socket
from Database import Database
import numpy as np


class Tss_map:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

        # self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293']
        self.cell_lines = ['K562', 'GM12878']
        self.chrom = ['chr1', 'chr2', 'chr3']

    def grouping_by_loc(self, df):
        df.loc[:, 'location'] = df['start'] + ';' + df['end']
        df_grp = df.groupby('location')

        df_res = []
        for group, df_sub in df_grp:
            df_res.append([*df_sub.iloc[0]['start': 'type'], ';'.join(df_sub['cell_line'])])
        return pd.DataFrame(df_res, columns=['start', 'end', 'type', 'cell_lines'])

    def fill_table(self, df, df_buffer):
        for idx in df.index:
            cell_lines = df.loc[idx, 'cell_lines'].split(';')
            location = '{};{}'.format(df.loc[idx, 'start'], df.loc[idx, 'end'])
            for cline in cell_lines:
                df_buffer.loc[location, cline] += 1
        return df_buffer

    def to_table(self):
        fpath = os.path.join(self.root, 'Papers/Tss_map', self.__class__.__name__ + '2.db')
        con = sqlite3.connect(fpath)

        out_path = os.path.join(self.root, 'Papers/Tss_map', self.__class__.__name__ + '.xlsx')
        N = len(self.cell_lines)
        dfs = []
        for chrom in self.chrom:
            print(chrom)
            for strand in ['+', '-']:
                tname = '{}_{}'.format(chrom, strand)
                df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
                index = df['start'].astype(str) + ';' + df['end'].astype(str)
                buffer = np.zeros((df.shape[0], N))

                df_buffer = pd.DataFrame(data=buffer, columns=self.cell_lines, index=index)
                df_buffer['type'] = df['type'].values
                dfs.append([self.fill_table(df, df_buffer), tname])

        with pd.ExcelWriter(out_path) as writer:
            for ele in dfs:
                df, tname = ele
                df.to_excel(writer, sheet_name=tname)

    def run(self):
        fpath = os.path.join(self.root, 'Papers/Tss_map', 'Tss_map_table' + '.db')
        con = sqlite3.connect(fpath)

        out_fpath = os.path.join(self.root, 'Papers/Tss_map', self.__class__.__name__ + '2.db')
        out_con = sqlite3.connect(out_fpath)

        for chrom in self.chrom:
            print(chrom)
            for strand in ['+', '-']:
                df_strand = []
                for cline in self.cell_lines:
                    tname = '{}_{}_{}'.format(cline, chrom, strand)
                    df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
                    df.loc[:, 'cell_line'] = cline
                    df_strand.append(df)

                df_strand = pd.concat(df_strand)
                df_res = self.grouping_by_loc(df_strand)
                df_res[['start', 'end']] = df_res[['start', 'end']].astype(int)
                df_res = df_res.sort_values(by='start')
                df_res.to_sql('{}_{}'.format(chrom, strand), out_con, if_exists='replace', index=None)


if __name__ == '__main__':
    tm = Tss_map()
    # tm.run()
    tm.to_table()
