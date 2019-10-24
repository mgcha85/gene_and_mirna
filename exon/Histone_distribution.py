import pandas as pd
import sqlite3
import os
import socket
from Database import Database
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Server import Server
import sys
import pickle as pkl
import numpy as np
import shutil
from exon.histogram_gpu import histogram_gpu


class Histone_distribution:
    def __init__(self, task_num=1, which='newton'):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.task_num = task_num
        self.which = which

    def run(self):
        out_path = os.path.join(self.root, 'database/RNA-seq/out/exon_analysis')
        flist = os.listdir(out_path)
        flist = [x for x in flist if x.endswith('.cha')]
        dfs = {}

        with open('intersection_genes.txt', 'rt') as f:
            genes = f.read().split('\n')

        hg = histogram_gpu()
        dbpath_src = ''
        src_tname = ''

        for fname in flist:
            fpath = os.path.join(out_path, fname)
            with open(fpath, 'rb') as f:
                dfs[os.path.splitext(fname)[0]] = pkl.load(f)
            print('complete loading {}'.format(fname))

            for gene in genes:
                df = dfs[fname]
                tissue = fname.replace('exon_distribution_', '')
                df_data = df.loc[gene, tissue]
                df_data = df_data[df_data['feature'] == 'transcript']
                df_data_grp = df_data.groupby('chromosome')
                dfs = {}
                for chr, df_data_sub in df_data_grp:
                    dfs[chr] = hg.histogram_gpu(df_data_sub, dbpath_src, src_tname)


