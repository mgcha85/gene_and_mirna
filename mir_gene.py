import sqlite3
import pandas as pd
import os
import socket
from scipy.stats import hypergeom
import numpy as np


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

    def test(self):
        fpath = os.path.join(self.root, 'database', 'important_genes_from_Amlan.xlsx')
        df_ref = pd.read_excel(fpath, index_col=0)
        df_ref = df_ref.dropna(subset=['genes'])

        df_high = pd.read_excel(os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_100.xlsx'))
        high_genes = set(df_high['gene_name'])

        for idx in df_ref.index:
            genes = set(df_ref.loc[idx, 'genes'].split(';'))
            compl = genes - high_genes
            print(len(compl))

    def comparison(self, hbw):
        fpath = os.path.join(self.root, 'database', 'important_genes_from_Amlan.db')
        con = sqlite3.connect(fpath)

        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_{}_2.xlsx'.format(hbw))
        writer = pd.ExcelWriter(out_path, engine='xlsxwriter')
        for tname in ['important_genes', 'mtb', 'ts']:
            df_ref = pd.read_sql("SELECT * FROM '{}' WHERE genes>''".format(tname), con, index_col='miRNA')
            genes = ';'.join(df_ref['genes'])
            genes = set(genes.split(';'))
            genes = [x for x in genes if len(x) > 0]
            m = len(genes)

            df_high = pd.read_excel(os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.xlsx'.format(hbw)))
            high_genes = set(df_high['gene_name'])
            n = len(high_genes) - m

            fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_{}.xlsx'.format(hbw))
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
                q = len(set.intersection(genes, genes_ref))
                k = len(genes)
                contents.append([mir, m, n, q, k])
            pd.DataFrame(contents, columns=['mir', 'm', 'n', 'q', 'k']).to_excel(writer, sheet_name=tname, index=None)
        writer.save()
        writer.close()

    def phypher(self, hbw):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_{}_2.xlsx'.format(hbw))

        writer = pd.ExcelWriter(fpath, engine='xlsxwriter')
        for tname in ['important_genes', 'mtb', 'ts']:
            df = pd.read_excel(fpath, index_col=0, sheet_name=tname)
            for idx in df.index:
                p = hypergeom.sf(df.loc[idx, 'q'], df.loc[idx, 'n'] + df.loc[idx, 'm'], df.loc[idx, 'm'], df.loc[idx, 'k'])
                df.loc[idx, 'phypher'] = -np.log(p)
            df.to_excel(writer, sheet_name=tname)
        writer.save()
        writer.close()

    def plot(self, hbw):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_{}_2.xlsx'.format(hbw))
        tnames = ['important_genes', 'mtb', 'ts']
        N = len(tnames)
        for i, tname in enumerate(tnames):
            df = pd.read_excel(fpath, index_col=0, sheet_name=tname)
            df['phypher'] = df['phypher'].replace(-np.inf, np.nan)
            df = df.dropna(subset=['phypher'])

            xaxis = range(df.shape[0])
            plt.subplot(N, 1, i+1)
            plt.scatter(xaxis, df['phypher'])
            plt.ylabel('log(phypher)')
            plt.xticks(xaxis, df.index, fontsize=6, rotation=30)
            plt.title('max: {:0.2f}, min: {:0.2f}, mean: {:0.2f}'.format(df['phypher'].max(), df['phypher'].min(), df['phypher'].mean()))
            # plt.savefig(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'log_phypher.png'))
        plt.show()


if __name__ == '__main__':
    mg = mir_gene()
    # mg.test()
    mg.comparison(100)
    mg.phypher(100)
    mg.plot(100)
