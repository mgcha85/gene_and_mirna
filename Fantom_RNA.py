import os
import sqlite3
import pandas as pd
import numpy as np
import socket
from Database import Database


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

    def get_tissues(self):
        df = pd.read_csv("E-MTAB-1733.csv")
        tissues = np.array([x.split('_')[0] for x in df['src']])
        tissues = np.unique(tissues)
        with open('tissues.txt', 'wt') as f:
            f.write('\n'.join(tissues))

    def extract_tissues_from_fantom(self):
        with open('tissues.txt', 'rt') as f:
            tissues = f.read().split('\n')

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df = pd.read_csv(fpath)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.tissues.db')
        con = sqlite3.connect(fpath)

        for tissue in tissues:
            print(tissue)
            columns = ['chromosome', 'start', 'end', 'strand', 'length']
            for col in df.columns[5:]:
                if tissue in col:
                    columns.append(col)
            df[columns].to_sql(tissue, con, if_exists='replace', index=None)

    def tpm(self, df_gene, df_rna):
        length = abs(df_gene['start'] - df_gene['end'])
        N = df_rna.shape[0]
        rpk = N / length
        return rpk / 1e6

    def rpkm(self, df_gene, df_rna):
        length = abs(df_gene['start'] - df_gene['end'])
        N = df_rna.shape[0]
        rpm = N / 1e6
        return rpm / length

    def merge_db(self):
        df = pd.read_csv("E-MTAB-1733.csv")
        new_col = [''] * df.shape[0]
        for i, url in enumerate(df.link):
            addr, fname = os.path.split(url)
            fid = os.path.splitext(fname)[0].split('_')[0]
            new_col[i] = fid
        df['fid'] = new_col
        df = df.drop_duplicates(subset=['fid'])
        df = df.set_index('fid', drop=True)

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq.db')
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_tissue.db'))

        tlist = Database.load_tableList(con)
        N = len(tlist)
        for i, fid in enumerate(tlist):
            print('{} / {}'.format(i + 1, N))
            tissue = df.loc[fid, 'src']

            df_sub = pd.read_sql_query("SELECT * FROM '{}'".format(fid), con)
            df_sub.to_sql(tissue, con_out, index=None, if_exists='append')

    def processInput(self, df, tissue, idx):
        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq.db')
        con_rna = sqlite3.connect(fpath)

        if idx % 1000 == 0 or idx + 1 == df.shape[0]:
            print('{:0.2f}%'.format(100 * (idx + 1) / df.shape[0]))

        chromosome = df.loc[idx, 'chromosome']
        start = df.loc[idx, 'start']
        end = df.loc[idx, 'end']
        strand = df.loc[idx, 'strand']

        df_rna = pd.read_sql_query("SELECT * FROM '{}_4a' WHERE chromosome='{}' AND strand='{}' AND NOT start>{end} "
                                   "AND NOT stop<{start}".format(tissue, chromosome, strand, start=start, end=end),
                                   con_rna)
        return self.rpkm(df.loc[idx], df_rna)

    def check_empty(self):
        df = pd.read_csv('E-MTAB-1733.csv')
        fid = []
        for idx in df.index:
            url = df.loc[idx, 'link']
            _, fname = os.path.split(url)
            fname, _ = os.path.splitext(fname)
            fname = fname.split('_')[0]
            fid.append(fname)
        fid = np.unique(np.array(fid))

        dirname = os.path.join(self.root, 'database/RNA-seq')
        flist = os.listdir(dirname)
        flist = [os.path.splitext(x)[0] for x in flist if x.endswith('.db')]

        remains = list(set(fid) - set(flist))
        print(remains)

    def run(self):
        from Database import Database
        from joblib import Parallel, delayed
        import multiprocessing
        # num_cores = 6
        num_cores = multiprocessing.cpu_count() // 2
        print('num_cores: {}'.format(num_cores))

        fpath_out = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.out.db')
        con_out = sqlite3.connect(fpath_out)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.tissues.db')
        con = sqlite3.connect(fpath)
        tissues = Database.load_tableList(con)

        for tissue in tissues:
            print(tissue)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con)
            if df.shape[1] <= 5:
                continue

            df['RPKM'] = Parallel(n_jobs=num_cores)(delayed(self.processInput)(df, tissue, idx) for idx in df.index)
            df.to_sql(tissue, con_out, index=None, if_exists='replace')
            break


if __name__ == '__main__':
    fr = Fantom_RNA()
    fr.run()
