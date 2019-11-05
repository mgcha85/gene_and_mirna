import sqlite3
import pandas as pd
import os
import socket
from scipy.stats import hypergeom


class mir_gene:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def to_sql(self):
        fpath = "/home/mingyu/Bioinformatics/database/feature_info.txt"
        df = pd.read_csv(fpath, sep='\t', names=['Target Gene', 'miRNA', 'binding energy', 'structural accessibility', 'context score'])
        df['Target Gene'] = df['Target Gene'].str.split('|', expand=True)[0]
        df['Target Gene'] = df['Target Gene'].str.split(' /// ', expand=True)[0]
        con = sqlite3.connect(fpath.replace('.txt', '.db'))
        df.to_sql('TargetScan', con, if_exists='replace', index=None)

    def TargetGene(self):
        fpath = os.path.join(self.root, 'database', 'feature_info.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'TargetScan'", con)
        df_grp = df.groupby('miRNA')

        contents = []
        for mir, df_sub in df_grp:
            tg = df_sub['Target Gene']
            contents.append([mir, ';'.join(tg)])
        df_res = pd.DataFrame(contents, columns=['miRNA', 'Target Gene'])
        df_res.to_excel('target_scan.xlsx', index=None)

    def mirTarbase(self):
        fpath = os.path.join(self.root, 'database', 'feature_info.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'miRTarBase' WHERE miRNA LIKE '%hsa%'", con)
        df_grp = df.groupby('miRNA')

        contents = []
        for mir, df_sub in df_grp:
            tg = df_sub['Target Gene']
            sm = df_sub['Species_miRNA']
            ex = df_sub['Experiments']
            contents.append([mir, ';'.join(tg), ';'.join(sm), ';'.join(ex)])
        df_res = pd.DataFrame(contents, columns=['miRNA', 'Target Gene', 'Species_miRNA', 'Experiments'])
        df_res.to_excel('miRTartBase.xlsx', index=None)

    def comparison(self):
        fpath = os.path.join(self.root, 'database', 'important_genes_from_Amlan.xlsx')
        df_ref = pd.read_excel(fpath, index_col=0)
        df_ref = df_ref.dropna(subset=['genes'])
        genes = ';'.join(df_ref['genes'])
        genes = set(genes.split(';'))

        df_high = pd.read_excel(os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_100.xlsx'))

        m = df_high.shape[0]
        common_genes = set.intersection(set(df_high['gene_name']) - genes)
        n = m - len(common_genes)

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_100.xlsx')
        df = pd.read_excel(fpath, index_col=0)

        mirs = set.intersection(set(df_ref.index), set(df.index))

        contents = []
        for mir in mirs:
            if isinstance(df.loc[mir, 'gene_name'], float):
                continue
            genes = set(df.loc[mir, 'gene_name'].split(';'))

            if isinstance(df_ref.loc[mir, 'genes'], float):
                continue
            genes_ref = set(df_ref.loc[mir, 'genes'].split(';'))
            contents.append([mir, m, n, len(set.intersection(genes, genes_ref)), len(genes)])
        pd.DataFrame(contents, columns=['mir', 'N', 'M', 'n', 'm']).to_excel(fpath.replace('.xlsx', '_2.xlsx'), index=None)

    def phyper(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_100_2.xlsx')
        df = pd.read_excel(fpath, index_col=0)

        for idx in df.index:
            p = hypergeom.sf(df.loc[idx, 'n'], df.loc[idx, 'N'] + df.loc[idx, 'M'], df.loc[idx, 'N'], df.loc[idx, 'm'])
            df.loc[idx, 'phyper'] = p
        df.to_excel(fpath)


if __name__ == '__main__':
    mg = mir_gene()
    mg.comparison()
    mg.phyper()
