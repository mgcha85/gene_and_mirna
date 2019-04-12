import pandas as pd
import sqlite3
import os
import socket
from Database import Database
import pickle as pkl


class Tss_map_table:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
            self.cell_lines = self.get_cell_lines()[300:]
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
            self.cell_lines = self.get_cell_lines()[300:]
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
            self.cell_lines = self.get_cell_lines()[:200]

        # self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293']
        self.chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                      'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

        # self.cell_lines = ['GM12878']
        # self.chrom = ['chr9', 'chr10', 'chr11', 'chr12']
        self.tnames = {'gencode': 'gencode_v28_transcripts_{}_{}'}

    def get_cell_lines(self):
        from Database import Database
        return Database.load_tableList(self.connect_cage_tags())

    def connect_gencode(self):
        fpath = os.path.join(self.root, 'database', 'genecode.db')
        return sqlite3.connect(fpath, check_same_thread=False)

    def connect_cage_tags(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc_full.db')
        return sqlite3.connect(fpath, check_same_thread=False)

    def connect_mirTSS(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        return sqlite3.connect(fpath, check_same_thread=False)

    def processInput_search_neighbor_tag(self, df_gencode, idx, tss_label, cline, strand):
        con_ctag = self.connect_cage_tags()
        tss = df_gencode.loc[idx, tss_label]
        gene = df_gencode.loc[idx, 'gene']
        chrom = df_gencode.loc[idx, 'chromosome']
        start = tss - 100
        end = tss + 100
        df_tags = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND strand='{}' AND NOT start>{end} AND "
                                    "NOT end<{start}".format(cline, chrom, strand, start=start, end=end), con_ctag)
        df_tags = df_tags[df_tags.iloc[:, 4:].sum(axis=1) > 0]
        if not df_tags.empty:
            df_tags = df_tags[['start', 'end']]
            if 'MIR' in gene:
                df_tags.loc[:, 'type'] = 'mTags'
                return df_tags
            else:
                df_tags.loc[:, 'type'] = 'gTags'
                return df_tags

    def search_neighbor_tag(self, df_gencode, tss_label, cline, strand, par=False):
        if par is False:
            tags = []
            for idx in df_gencode.index:
                tags.append(self.processInput_search_neighbor_tag(df_gencode, idx, tss_label, cline, strand))
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            num_cores = multiprocessing.cpu_count()
            if num_cores > 10:
                num_cores = 10

            tags = Parallel(n_jobs=num_cores)(delayed(self.processInput_search_neighbor_tag)(df_gencode, idx, tss_label, cline, strand) for idx in df_gencode.index)
        return pd.concat(tags)

    def get_others(self, df, cline, chrom):
        con_ctag = self.connect_cage_tags()
        df.index = df['start'].astype(str) + ';' + df['end'].astype(str)

        df_tags = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}'".format(cline, chrom), con_ctag)
        df_tags.index = df_tags['start'].astype(str) + ';' + df_tags['end'].astype(str)

        others = set(df_tags.index) - set(df.index)
        df_tags = df_tags.loc[others]
        df_tags.loc[:, 'type'] = 'oTag'
        df = pd.concat([df, df_tags[df.columns]])
        return df.sort_values(by='start')

    def set_type(self, df):
        df['location'] = df['start'].astype(str) + ';' + df['end'].astype(str)
        df_grp = df.groupby('location')

        df_res = []
        for group, df_sub in df_grp:
            df_sub = df_sub.drop_duplicates()
            if df_sub.shape[0] > 1:
                df_res.append([*df_sub.iloc[0][['start', 'end']].values.flatten(), ';'.join(df_sub['type'])])
            else:
                df_res.append([*df_sub[['start', 'end', 'type']].values.flatten()])
        return pd.DataFrame(df_res, columns=['start', 'end', 'type'])
        
    def search_table(self):
        con_gencode = self.connect_gencode()
        # con_mtss = self.connect_mirTSS()
        # con_ctag = self.connect_cage_tags()

        fpath_out = os.path.join(self.root, 'Papers/Tss_map', self.__class__.__name__ + '.db')
        print(fpath_out)
        con_out = sqlite3.connect(fpath_out, check_same_thread=False)

        for i, cline in enumerate(self.cell_lines):
            print(i, cline)
            for chrom in self.chrom:
                for slabel, strand, tss_label in zip(['plus', 'minus'], ['+', '-'], ['start', 'end']):
                    print(chrom, strand)
                    df_gencode = pd.read_sql_query("SELECT * FROM '{}'".format(self.tnames['gencode'].format(chrom, slabel)), con_gencode)
                    df_ctag = self.search_neighbor_tag(df_gencode, tss_label, cline, strand)

                    df_ctag = self.set_type(df_ctag)
                    df_ctag = self.get_others(df_ctag, cline, chrom)
                    df_ctag.to_sql('{}_{}_{}'.format(cline, chrom, strand), con_out, index=None, if_exists='replace')


if __name__ == '__main__':
    tm = Tss_map_table()
    tm.search_table()

