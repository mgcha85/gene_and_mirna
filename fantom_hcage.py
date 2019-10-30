import os
import sqlite3
import pandas as pd
from Database import Database
import math


class StdevFunc:
    def __init__(self):
        self.M = 0.0
        self.S = 0.0
        self.k = 1

    def step(self, value):
        if value is None:
            return
        tM = self.M
        self.M += (value - tM) / self.k
        self.S += (value - tM) * (value - self.M)
        self.k += 1

    def finalize(self):
        if self.k < 3:
            return None
        return math.sqrt(self.S / (self.k-2))


class Fantom:
    def __init__(self):
        import socket
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.cells = []
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    def check_score(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'FANTOM_tissue.db')
        con = sqlite3.connect(fpath)
        con.create_aggregate("STDEV", 1, StdevFunc)

        tlist = Database.load_tableList(con)
        dcols = ['chromosome', 'start', 'end', 'score', 'strand']

        df_res = pd.DataFrame(index=tlist, columns=['SUM', 'AVG', 'STDEV'])
        for tname in tlist:
            print(tname)
            columns = Database.load_table_columns(con, tname)
            for d in dcols:
                columns.remove(d)

            for cal in ['SUM', 'AVG', 'STDEV']:
                cols = []
                for c in columns:
                    cols.append('{}("{}")'.format(cal, c))

                sql = "SELECT {} FROM '{}'".format(', '.join(cols), tname)
                df = pd.read_sql(sql, con)
                stats = df.iloc[0, :].astype(int)
                df_res.loc[tname, cal] = ','.join(stats.astype(str))
        df_res.to_excel(os.path.join(self.root, 'database/Fantom/v5/tissues', 'FANTOM_check_score.xlsx'))

    def check_score_celllines(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        fpath = os.path.join(dirname, 'human_hCAGE_celllines.db')
        con = sqlite3.connect(fpath)

        tlist = Database.load_tableList(con)
        df = pd.Series()
        for tname in tlist:
            print(tname)
            score = pd.read_sql("SELECT SUM(score) FROM '{}'".format(tname), con)
            df[tname] = score.loc[0, 'SUM(score)']
        df.to_csv(os.path.join(dirname, 'human_hCAGE_celllines.csv'))

    def subset_test(self):
        dirname = os.path.join(self.root, 'database/gencode')
        dfs = {}
        for hbw in [100, 300, 500]:
            fpath = os.path.join(dirname, 'high_correlated_fan_rna_{}.xlsx'.format(hbw))
            dfs[hbw] = pd.read_excel(fpath, index_col=7)

        from itertools import combinations
        for comb in combinations([100, 300, 500], 2):
            set1 = set(dfs[comb[0]].index)
            set2 = set(dfs[comb[1]].index)
            inters = set.intersection(set1, set2)
            comp_1 = set1 - inters
            comp_2 = set2 - inters

            print(comb, dfs[comb[0]].shape[0], dfs[comb[1]].shape[0], len(inters), len(comp_1), len(comp_2))


if __name__ == '__main__':
    f = Fantom()
    # f.check_score()
    f.check_score_celllines()
