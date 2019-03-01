import pandas as pd
import sqlite3
import os
import socket


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

        self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293']
        self.chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                      'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM',
                      'chrX', 'chrY']
        self.tnames = {'gencode': 'gencode_v28_transcripts_{}_{}'}

    def connect_gencode(self):
        fpath = os.path.join(self.root, 'database', 'genecode.db')
        return sqlite3.connect(fpath)

    def connect_cage_tags(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc.db')
        return sqlite3.connect(fpath)

    def connect_mirTSS(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        return sqlite3.connect(fpath)

    def processInput_search_neighbor_tag(self, df_gencode, idx, tss_label, cline):
        con_ctag = self.connect_cage_tags()
        tss = df_gencode.loc[idx, tss_label]
        start = tss - 100
        end = tss + 100
        df_tags = pd.read_sql_query("SELECT * FROM {} WHERE NOT start>{end} AND NOT end<{start}"
                                    "".format(cline, start=start, end=end), con_ctag)
        df_tags = df_tags[df_tags.iloc[:, 4:].sum(axis=1) > 0]
        if not df_tags.empty:
            return df_tags[['start', 'end']]

    def search_neighbor_tag(self, df_gencode, tss_label, cline):
        from joblib import Parallel, delayed
        import multiprocessing

        num_cores = multiprocessing.cpu_count()
        tags = Parallel(n_jobs=num_cores)(delayed(self.processInput_search_neighbor_tag)(df_gencode, idx, tss_label, cline) for idx in df_gencode.index)

        # for idx in df_gencode.index:
        #     tss = df_gencode.loc[idx, tss_label]
        #     start = tss - 100
        #     end = tss + 100
        #     df_tags = pd.read_sql_query("SELECT * FROM {} WHERE NOT start>{end} AND NOT end<{start}"
        #                                 "".format(cline, start=start, end=end), con_ctag)
        #     df_tags = df_tags[df_tags.iloc[:, 4:].sum(axis=1) > 0]
        #     if not df_tags.empty:
        #         tags.append(df_tags[['start', 'end']])
        return pd.concat(tags)

    def run(self):
        con_gencode = self.connect_gencode()
        con_mtss = self.connect_mirTSS()
        # con_ctag = self.connect_cage_tags()

        fpath_out = os.path.join(self.root, 'Papers/Tss_map', self.__class__.__name__ + '.db')
        print(fpath_out)
        con_out = sqlite3.connect(fpath_out)

        for cline in self.cell_lines:
            for chrom in self.chrom:
                for slabel, strand, tss_label in zip(['plus', 'minus'], ['+', '-'], ['start', 'end']):
                    print(cline, chrom, strand)
                    df_gencode = pd.read_sql_query("SELECT start, end FROM '{}'".format(self.tnames['gencode'].format(chrom, slabel)), con_gencode)
                    df_mtss = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE chromosome='{}' "
                                                "AND strand='{}'".format(chrom, strand), con_mtss)

                    df_ctag = self.search_neighbor_tag(df_gencode, tss_label, cline)

                    df_mtss = df_mtss[['tss', 'tss']]
                    df_mtss.columns = ['start', 'end']

                    df_ctag.loc[:, 'type'] = 'tag'
                    df_mtss.loc[:, 'type'] = 'miRNA'
                    df_gencode.loc[:, 'type'] = 'gene'

                    df_res = pd.concat([df_ctag, df_mtss, df_gencode])
                    df_res.to_sql('{}_{}_{}'.format(cline, chrom, strand), con_out, index=None)


if __name__ == '__main__':
    tm = Tss_map()
    tm.run()
