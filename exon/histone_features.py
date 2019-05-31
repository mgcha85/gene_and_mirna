import os
import sqlite3
import pandas as pd
import socket
import pickle as pkl

from XmlHandler import XmlHandler
from histogram_gpu import histogram_gpu


class histone_features:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
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
                df = hgpu.histogram_gpu(df_ref, rsc_path, fname)
                dfs[histone] = df
                return

            # for idx in df_ref.index:
            #     gene_name = df_ref.loc[idx, 'gene_name']
            #     print(gene_name)



if __name__ == '__main__':
    hf = histone_features()
    hf.run()