import os
import sqlite3
import pandas as pd
import socket
import pickle as pkl
import numpy as np
from XmlHandler import XmlHandler
from histogram_gpu import histogram_gpu


class histone_features:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.histones = {'H3K4me1': 0, 'H3K4me2': 1, 'H3K4me3': 2, 'H3K9ac': 3, 'H3K9me3': 4, 'H3K27ac': 5,
                         'H3K27me3': 6, 'H3K36me3': 7}

    def get_rna_list(self):
        dirname = os.path.join(self.root, 'database/RNA-seq/out/exon_analysis')
        flist = os.listdir(dirname)

        fid = []
        for fname in flist:
            fname = fname.replace('exon_distribution_', '')
            fid.append(os.path.splitext(fname)[0])
        return fid

    def load_ref(self, tissue):
        fid = self.get_rna_list()
        fid = [x for x in fid if tissue in x]
        if len(fid) == 0:
            return

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql_query("SELECT * FROM '{}' WHERE gene_name>'' AND NOT chromosome='chrMT' AND "
                                 "NOT chromosome='chrY'".format(fid[0]), con)
        return df[df['chromosome'].str.len() <= 5].reset_index(drop=True)

    def get_hist(self, df_ref, rsc_path, fname, BINSIZE=100):
        con = sqlite3.connect(rsc_path)

        hist_out = []
        for idx in df_ref.index:
            chromosome = df_ref.loc[idx, 'chromosome']
            start = df_ref.loc[idx, 'start']
            end = df_ref.loc[idx, 'end']
            df__ = pd.read_sql_query("SELECT * FROM '{}_{}' WHERE NOT start>{end} AND NOT end<{start}".format(fname, chromosome, start=start, end=end), con)

            bin_width = (end - start + BINSIZE) // BINSIZE
            hist = np.zeros(bin_width)
            if df__.empty:
                hist_out.append(hist)
                continue

            bin_start = start
            for i in range(bin_width):
                bin_end = bin_start + BINSIZE
                df = df__[~(~(df__['start'] > bin_end) & ~(df__['end'] < bin_start))]
                hist[i] = df.shape[0]
                bin_start += BINSIZE
            hist_out.append(hist)
        return hist_out

    def run(self):
        import matplotlib.pyplot as plt

        dirname__ = os.path.join(self.root, 'database/Histone ChIP-seq')
        filenames = XmlHandler.load_param('histogram_param.xml')
        N = len(self.histones)
        hgpu = histogram_gpu(fpath='histogram_param.xml')

        for tissue, fnames in filenames.items():
            df_ref = self.load_ref(tissue)
            if df_ref is None:
                continue

            dfs = {}
            for histone, fname in fnames.items():
                rsc_path = os.path.join(dirname__, tissue, 'histone.db')
                hist_out = self.get_hist(df_ref, rsc_path, fname)
                print(hist_out)

            # for idx in df_ref.index:
            #     gene_name = df_ref.loc[idx, 'gene_name']
            #     print(gene_name)



if __name__ == '__main__':
    hf = histone_features()
    hf.run()