import socket
import os
import sqlite3
import pandas as pd
import numpy as np
import shutil


class compare_conv_files:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinfomatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def load_ref(self):
        with open('temp.csv', 'rt') as f:
            contents = f.read().split('\n')

        fid = np.array([None]*len(contents))
        for i, line in enumerate(contents):
            ulr, fname_ = os.path.split(line)
            fname = fname_.split('.')[0]
            fid[i] = fname.split('_')[0]
        return np.unique(fid)

    def get_remain(self):
        fid = self.load_ref()

        dirname = os.path.join(self.root, 'database/RNA-seq/gtf')
        flist = os.listdir(dirname)

        fnames = []
        for fname in flist:
            fname, ext = os.path.splitext(fname)
            fnames.append(fname)

        return sorted(list(set(fid) - set(fnames)))

    def move_fastaq(self):
        remain = self.get_remain()

        dirname = os.path.join(self.root, 'database/RNA-seq/fastq')
        flist = [x for x in os.listdir(dirname) if x.endswith('.gz')]
        fpaths = {}
        for fname in flist:
            fid = fname.split('.')[0]
            fid = fid.split('_')[0]

            if fid in fpaths:
                fpaths[fid].append(os.path.join(dirname, fname))
            else:
                fpaths[fid] = [os.path.join(dirname, fname)]

        N = len(remain)
        M = (N + 8) // 9
        for i, fid in enumerate(remain):
            dirnum = (i // M) + 1
            for fpath in fpaths[fid]:
                dirname, fname = os.path.split(fpath)
                shutil.move(fpath, os.path.join(dirname, str(dirnum), fname))


if __name__ == '__main__':
    ccf = compare_conv_files()
    print(ccf.get_remain())
