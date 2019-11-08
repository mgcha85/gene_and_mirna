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
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.table_names = {}

    def to_server(self):
        from Server import Server
        import sys

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
        sql = "SELECT chromosome, start, end, strand, FPKM, reference_id, ref_gene_id, ref_gene_name FROM '{}'".format(tname)
        # sql = "SELECT chromosome, start, end, strand, score FROM '{}'".format(tname)
        df = pd.read_sql_query(sql, con)

        out_con = sqlite3.connect(out_path)
        for chr, df_chr in df.groupby('chromosome'):
            if len(chr) > 5:
                continue
            for str, df_str in df_chr.groupby('strand'):
                df_str = df_str.sort_values(by=['start', 'end'])
                df_str.drop(['chromosome', 'strand'], axis=1).to_sql('_'.join([tname, chr, str]), out_con, if_exists='replace', index=None)


if __name__ == '__main__':
    ut = Util()
    if ut.hostname == 'mingyu-Precision-Tower-7810':
        ut.to_server()
    else:
        # fpath = os.path.join(ut.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        fpath = os.path.join(ut.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        tlist = Database.load_tableList(sqlite3.connect(fpath))
        for tname in tlist:
            ut.split(fpath, tname)
