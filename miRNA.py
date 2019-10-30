import pandas as pd
import numpy as np
import os
import sqlite3
import socket


class miRNA:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def convert(self):
        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.bed")
        df = pd.read_csv(fpath, sep='\t', index_col=0, names=['locations'])
        df.index.name = 'miRNA'

        out_path = fpath.replace('.bed', '.db')
        con = sqlite3.connect(out_path)
        df.to_sql('original', con)

    def seperate_loc(self):
        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.db")
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'original'", con)
        df['locations'] = df['locations'].str.split(';')

        df_res = pd.DataFrame(index=df.index, columns=['miRNA', 'chromosome', 'starts', 'ends'])
        for idx in df.index:
            mir = df.loc[idx, 'miRNA']
            loc = df.loc[idx, 'locations']
            N = len(loc)
            df_each = pd.DataFrame(data=np.zeros((N, 3)), columns=['chromosome', 'start', 'end'], index=range(N))

            for i, l in enumerate(loc):
                chr, start, end = l.split('_')
                df_each.loc[i, 'chromosome'] = chr
                df_each.loc[i, 'start'] = int(start)
                df_each.loc[i, 'end'] = int(end)

            df_res.loc[idx, 'miRNA'] = mir
            df_res.loc[idx, 'chromosome'] = chr
            df_res.loc[idx, 'starts'] = ';'.join(df_each['start'].astype(int).astype(str))
            df_res.loc[idx, 'ends'] = ';'.join(df_each['end'].astype(int).astype(str))
        df_res.to_sql('seperate_tss', con, index=None)

    def mean_tss(self):
        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.db")
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'seperate_tss'", con, index_col='miRNA')
        df = df.loc[~df.index.duplicated(keep='first')]

        df_res = pd.DataFrame(index=df.index, columns=['chromosome', 'tss'])
        for idx in df.index:
            starts = df.loc[idx, 'starts'].split(';')
            ends = df.loc[idx, 'ends'].split(';')

            smin = min(map(int, starts))
            emin = min(map(int, ends))
            smax = max(map(int, starts))
            emax = max(map(int, ends))

            tss = (min([smin, emin]) + max([smax, emax])) // 2
            df_res.loc[idx, 'chromosome'] = df.loc[idx, 'chromosome']
            df_res.loc[idx, 'tss'] = tss
        df_res.to_sql('mean_tss', con)

    def add_strand(self):
        from Database import Database

        con_ref = sqlite3.connect(os.path.join(self.root, 'database', 'papers.db'))
        tlist = Database.load_tableList(con_ref)
        tlist = [x for x in tlist if x != 'miRNA_Type']

        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.db")
        con = sqlite3.connect(fpath)

        df = pd.read_sql("SELECT * FROM 'mean_tss'", con)
        for idx in df.index:
            mirname = df.loc[idx, 'miRNA']
            print(mirname)

            for tname in tlist:
                columns = Database.load_table_columns(con_ref, tname)
                try:
                    sidx = columns.index('strand')
                except:
                    sidx = columns.index('Strand')

                df_ref = pd.read_sql("SELECT miRNA, {} FROM '{}' WHERE miRNA='{}'".format(columns[sidx], tname, mirname), con_ref, index_col='miRNA')
                if not df_ref.empty:
                    if df_ref.shape[0] > 1:
                        df_ref = df_ref.iloc[:1]
                    strand = df_ref.loc[mirname, columns[sidx]]
                    df.loc[idx, 'strand'] = strand
                    break
        df.to_sql('mean_tss', con, if_exists='replace', index=None)


if __name__ == '__main__':
    mir = miRNA()
    # mir.convert()
    # mir.seperate_loc()
    # mir.mean_tss()
    mir.add_strand()
