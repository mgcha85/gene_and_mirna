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


class Exon_distribution:
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

    def to_server(self):
        server = Server(self.root)
        which = 'newton'
        server.connect(which=which)

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='04:00:00', which=which)
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def get_intersection_genes(self):
        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_full.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        genes = []
        for tname in tlist:
            print(tname)
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            if df.empty:
                continue
            genes.append(set(df['gene_name']))

        with open('temp.cha', 'wb') as f:
            pkl.dump(genes, f)

        genes = list(set.intersection(*genes))
        genes = [x for x in genes if x is not None]
        with open('intersection_genes.txt', 'wt') as f:
            f.write('\n'.join(genes))

    def select_one_trans(self, df):
        df_tr = df[df['feature'] == 'transcript']
        if df_tr.shape[0] == 1:
            return df
        else:
            df_tr = df_tr.reset_index()
            midx = df_tr['cov'].astype(float).idxmax()
            sidx = df_tr.loc[midx, 'index']
            if midx + 1 >= df_tr.shape[0]:
                eidx = df_tr.shape[0]
            else:
                eidx = df_tr.loc[midx+1, 'index']
            return df.loc[sidx: eidx - 1, :]

    def choose_one_rep(self):
        dirname = os.path.join(self.root, 'database/RNA-seq/out/exon_analysis')
        flist = os.listdir(dirname)
        flist = [x for x in flist if x.endswith('.cha')]

        contents = {}
        for fname__ in flist:
            fname = os.path.splitext(fname__)[0]
            tissue = fname.replace('exon_distribution_', '')
            tissue = tissue.split('_')[0]
            if tissue not in contents:
                contents[tissue] = [fname__]
            else:
                contents[tissue].append(fname__)

        for tissue, fnames in contents.items():
            for i, fname in enumerate(fnames):
                super_dir = '/'.join(dirname.split('/')[:-1])
                if i > 0:
                    shutil.move(os.path.join(dirname, fname), os.path.join(super_dir, fname))

    def plot(self, gene, flist, dfs):
        print(gene)
        plt.figure(figsize=(12, 8))
        for j, fname in enumerate(flist):
            fname = os.path.splitext(fname)[0]
            df = dfs[fname]
            tissue = fname.replace('exon_distribution_', '')
            df_data = df.loc[gene, tissue]
            df_data = self.select_one_trans(df_data)
            for sidx in df_data.index:
                feature = df_data.loc[sidx, 'feature']
                if feature == 'transcript':
                    color = 'g'
                    ylen = 1
                else:
                    color = 'c'
                    ylen = 0.5

                chromosome = df_data.loc[sidx, 'chromosome']
                start = df_data.loc[sidx, 'start']
                end = df_data.loc[sidx, 'end']
                strand = df_data.loc[sidx, 'strand']
                cov = df_data.loc[sidx, 'cov']
                length = end - start

                plt.subplot(5, 6, j + 1)
                currentAxis = plt.gca()
                currentAxis.add_patch(Rectangle((start, 0), length, ylen, alpha=0.5, color=color))
                if feature == 'transcript':
                    xaxis = np.linspace(start, end, 4).astype(int)
                    plt.xticks(xaxis, xaxis.astype(str), rotation=30, fontsize=6)
                    plt.xlim([start - 500, end + 500])
                    ttl = ':'.join([chromosome, str(start), str(end), strand, '{:0.4f}'.format(float(cov))])

            if j == 0:
                plt.title(tissue + ' [' + ttl + ']', fontsize=6)
            else:
                plt.title(tissue, fontsize=6)
        out_path = os.path.join(self.root, 'database/RNA-seq/out/figures', gene + '.png')
        plt.savefig(out_path)
        plt.close()

    def get_distribution(self):
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores > 6:
            num_cores = 6

        out_path = os.path.join(self.root, 'database/RNA-seq/out/exon_analysis')
        flist = os.listdir(out_path)
        flist = [x for x in flist if x.endswith('.cha')]
        dfs = {}

        for fname in flist:
            fpath = os.path.join(out_path, fname)
            with open(fpath, 'rb') as f:
                dfs[os.path.splitext(fname)[0]] = pkl.load(f)
            print('complete loading {}'.format(fname))

        with open('intersection_genes.txt', 'rt') as f:
            genes = f.read().split('\n')

        Parallel(n_jobs=num_cores)(delayed(self.plot)(gene, flist, dfs) for gene in genes)
        # for i, gene in enumerate(genes):
        #     print(gene)
        #
        #     plt.figure(figsize=(12, 8))
        #     for j, fname in enumerate(flist):
        #         fname = os.path.splitext(fname)[0]
        #         df = dfs[fname]
        #         tissue = fname.replace('exon_distribution_', '')
        #         df_data = df.loc[gene, tissue]
        #         df_data = self.select_one_trans(df_data)
        #         for sidx in df_data.index:
        #             feature = df_data.loc[sidx, 'feature']
        #             if feature == 'transcript':
        #                 color = 'g'
        #                 ylen = 1
        #             else:
        #                 color = 'c'
        #                 ylen = 0.5
        #
        #             chromosome = df_data.loc[sidx, 'chromosome']
        #             start = df_data.loc[sidx, 'start']
        #             end = df_data.loc[sidx, 'end']
        #             strand = df_data.loc[sidx, 'strand']
        #             cov = df_data.loc[sidx, 'cov']
        #             length = end - start
        #
        #             plt.subplot(4, 6, j + 1)
        #             currentAxis = plt.gca()
        #             currentAxis.add_patch(Rectangle((start, 0), length, ylen, alpha=0.5, color=color))
        #             if feature == 'transcript':
        #                 xaxis = np.linspace(start, end, 4).astype(int)
        #                 plt.xticks(xaxis, xaxis.astype(str), rotation=30, fontsize=6)
        #                 plt.xlim([start - 500, end + 500])
        #                 ttl = ':'.join([chromosome, str(start), str(end), strand, '{:0.4f}'.format(float(cov))])
        #
        #         if j == 0:
        #             plt.title(tissue + ' [' + ttl + ']', fontsize=6)
        #         else:
        #             plt.title(tissue, fontsize=6)
        #     out_path = os.path.join(self.root, 'database/RNA-seq/out/figures', gene + '.png')
        #     plt.savefig(out_path)
        #     plt.close()

    def process_run(self, tname, gene, i, N):
        if i % 1000 == 0 or i + 1 == N:
            print('{:0.2f}%'.format(100 * (i + 1) / N))
        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_full.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM '{}' WHERE gene_name='{}'".format(tname, gene), con)
        if df.empty:
            return
        return df

    def divide_tasks(self, num_tasks):
        with open('tissues__.txt', 'rt') as f:
            tlist__ = f.read().split('\n')

        N = len(tlist__)
        tlist = []
        M = (N + num_tasks - 1) // num_tasks
        for i in range(num_tasks):
            offset = i * M
            if i + M >= N:
                tlist.append(tlist__[offset:])
            else:
                tlist.append(tlist__[offset: offset + M])

        for i in range(num_tasks):
            tlist_batch = tlist[i]
            with open('{}/tissues.txt'.format(i + 1), 'wt') as f:
                f.write('\n'.join(tlist_batch))

    def run(self):
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores > 6:
            num_cores = 6

        with open('intersection_genes.txt', 'rt') as f:
            genes = f.read().split('\n')

        with open('tissues.txt', 'rt') as f:
            tlist = f.read().split('\n')

        print(tlist)
        N = len(genes)
        for tname in tlist:
            print(tname)
            contents = Parallel(n_jobs=num_cores)(delayed(self.process_run)(tname, gene, i, N) for i, gene in enumerate(genes))

            out_path = os.path.join(self.root, 'database/RNA-seq/out', 'exon_distribution_{}.cha'.format(tname))
            df = pd.DataFrame(contents, index=genes, columns=[tname])
            with open(out_path, 'wb') as f:
                pkl.dump(df, f)


if __name__ == '__main__':
    ed = Exon_distribution()
    if ed.hostname == 'mingyu-Precision-Tower-7810':
        ed.to_server()
    else:
        ed.get_distribution()
