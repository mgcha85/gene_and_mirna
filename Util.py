import sqlite3
import pandas as pd
from Database import Database
from Server import Server
import sys
import os
import socket


class Util:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.table_names = {}

    def to_server(self):
        which = 'newton'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='04:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def split(self, fpath, tname):
        out_path = fpath.replace('.db', '_spt.db')

        con = sqlite3.connect(fpath)
        df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)

        out_con = sqlite3.connect(out_path)
        for chr, df_chr in df.groupby('chromosome'):
            for str, df_str in df_chr.groupby('strand'):
                df_str = df_str.sort_values(by=['start', 'end'])
                df_str.to_sql('_'.join([tname, chr, str]), out_con, if_exists='replace', index=None)


if __name__ == '__main__':
    ut = Util()
    if ut.hostname == 'mingyu-Precision-Tower-7810':
        ut.to_server()
    else:
        # fpath = os.path.join(ut.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        fpath = os.path.join(ut.root, 'database/Fantom/v5/tissues', 'fantom_cage_by_tissue.db')
        tlist = Database.load_tableList(sqlite3.connect(fpath))
        for tname in tlist:
            ut.split(fpath, tname)
