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


class Exon_distribution:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def to_server(self):
        server = Server()
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server_root = os.path.join(server.server, 'source/gene_and_mirna/exon')
        server.job_script(fname, src_root=server_root, time='06:00:00')

        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command(
            "cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def get_intersection_genes(self):
        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_full.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        genes = []
        for tname in tlist:
            print(tname)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            if df.empty:
                continue
            genes.append(set(df['gene_name']))

        with open('temp.cha', 'wb') as f:
            pkl.dump(genes, f)

        genes = list(set.intersection(*genes))
        genes = [x for x in genes if x is not None]
        with open('intersection_genes.txt', 'wt') as f:
            f.write('\n'.join(genes))

    def get_distribution(self):
        out_path = os.path.join(self.root, 'database/RNA-seq/out', 'exon_distribution.cha')
        with open(out_path, 'rb') as f:
            df = pkl.load(f)

        for idx in df.index:
            plt.figure(figsize=(10, 8))

            for i, col in enumerate(df.columns):
                df_data = df.loc[idx, col]
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

                    plt.subplot(4, 5, i + 1)
                    currentAxis = plt.gca()
                    currentAxis.add_patch(Rectangle((start, 0), length, ylen, alpha=0.5, color=color))
                    if feature == 'transcript':
                        xaxis = np.linspace(start, end, 4).astype(int)
                        plt.xticks(xaxis, xaxis.astype(str), rotation=30, fontsize=6)
                        plt.xlim([start - 500, end + 500])
                        ttl = ':'.join([chromosome, str(start), str(end), strand, '{:0.4f}'.format(float(cov))])
                        plt.title(ttl, fontsize=6)

            out_path = os.path.join(self.root, 'database/RNA-seq/out/figures', idx + '.png')
            plt.savefig(out_path)

    def process_run(self, tname, gene, i, N):
        if i % 1000 == 0 or i + 1 == N:
            print('{:0.2f}%'.format(100 * (i + 1) / N))
        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_full.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql_query("SELECT * FROM '{}' WHERE gene_name='{}'".format(tname, gene), con)
        if df.empty:
            return
        return df

    def run(self):
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores > 6:
            num_cores = 6

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_full.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        with open('intersection_genes.txt', 'rt') as f:
            genes = f.read().split('\n')

        contents = []
        N = len(genes)
        for tname in tlist:
            print(tname)
            contents.append(Parallel(n_jobs=num_cores)(
                delayed(self.process_run)(tname, gene, i, N) for i, gene in enumerate(genes)))

        out_path = os.path.join(self.root, 'database/RNA-seq/out', 'exon_distribution.cha')
        df = pd.DataFrame(contents, index=tlist, columns=genes).T
        with open(out_path, 'wb') as f:
            pkl.dump(df, f)


if __name__ == '__main__':
    ed = Exon_distribution()
    if ed.hostname == 'mingyu-Precision-Tower-7810':
        # ed.get_intersection_genes()
        # ed.run()
        # ed.get_distribution()
        ed.to_server()
    else:
        ed.get_intersection_genes()
        ed.run()
        ed.get_distribution()
