import sqlite3
import pandas as pd
import os
import socket


class answer:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def union_sinigificant_mir(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation/nz')
        fpath = os.path.join(dirname, 'cross_stats.xlsx')
        reader = pd.ExcelFile(fpath)

        mir = []
        for sname in reader.sheet_names:
            df = reader.parse(sname)
            mir.append(set(df[df['median_diff_ratio'] < 1.1]['miRNA (test)']))

        mir = set.union(*mir)
        with open(os.path.join(dirname, 'significant_mir_union.txt'), 'wt') as f:
            f.write('\n'.join(mir))
        print(len(mir))

    def identical_mir(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation/nz')
        fpath = os.path.join(dirname, 'cross_stats.xlsx')
        reader = pd.ExcelFile(fpath)

        df_tis = pd.read_excel(os.path.join(self.root, 'database/Fantom/v5/tissues/out/cross_validation/nz', 'cross_stats.xlsx'), index_col=0)

        dfs = []
        for sname in reader.sheet_names:
            df = reader.parse(sname, index_col=0)
            dfs.append(df['median_diff_ratio'])
        df_res = pd.concat(dfs, axis=1)
        df_res.columns = range(10)
        df_res['mean'] = df_res.mean(axis=1)
        df_res['std'] = df_res.std(axis=1)
        df_res.loc[df_tis.index, 'tissue'] = df_tis['median_diff_ratio']
        df_res.sort_values(by='std').to_excel(os.path.join(dirname, 'identical_mir.xlsx'))

    def compare_go(self):
        with open(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation/nz', 'significant_mir_union.txt')) as f:
            sig_mir = set(f.read().split('\n'))

        df_go = pd.read_excel(os.path.join(self.root, 'database/target_genes', 'hyper_test_lasso_100_nz.xlsx'), index_col=0)
        df_go = df_go[(df_go['# significant'] > 0) | (df_go['q significant'] < 1)]
        mir = set(df_go.index)

        inter = set.intersection(sig_mir, mir)
        union = set.union(sig_mir, mir)
        with open(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation/nz', 'compare_go.txt'), 'wt') as f:
            f.write('\n'.join(inter))
        print(len(inter), len(union))


if __name__ == '__main__':
    ans = answer()
    ans.compare_go()
