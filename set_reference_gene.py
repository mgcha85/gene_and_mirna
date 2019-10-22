import pandas as pd
import os
import sqlite3
import socket
import csv


class set_reference_gene:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinfomatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.cells = []
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    def run(self):
        fpath = os.path.join(self.root, 'software/stringtie-1.3.3b', 'genes.gtf.gz')
        df = pd.read_csv(fpath, compression='gzip', names=self.gtf_columns, comment='#', sep='\t')
        df = df[~((df['feature'] == 'transcript') & (df['source'] == 'transcribed_unprocessed_pseudogene'))]
        df = df[df['source'] != 'transcribed_unprocessed_pseudogene']
        df.to_csv(fpath.replace('.gz', ''), sep='\t', index=None, header=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    srg = set_reference_gene()
    srg.run()
