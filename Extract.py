import sqlite3
import os
import pandas as pd
import numpy as np
import socket
import sys
import multiprocessing
from Server import Server


class Extract:
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
        self.num_cores = multiprocessing.cpu_count()

    def to_server(self):
        server = Server()
        which = 'newton'
        server.connect(which=which)

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server.job_script(fname, time='04:00:00', which=which)

        server_root = os.path.join(server.server, 'source/gene_and_mirna')
        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def get_high_corr(self, df):
        contents = []
        for mir in df.index:
            row = df.loc[mir, :]
            row = row[row > 0.6]
            if len(row) > 3:
                contents.append([mir, ';'.join(row.index), ';'.join(row.round(2).astype(str))])
        return pd.DataFrame(data=contents, columns=['miRNA', 'GENEs', 'corr (pearson)'])

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue_100_score_vector_corr.xlsx')
        df = pd.read_excel(fpath, index_col=0).T
        df = self.get_high_corr(df)
        df.to_excel(fpath.replace('.xlsx', '2.xlsx'), index=None)


if __name__ == '__main__':
    ex = Extract()
    ex.run()
