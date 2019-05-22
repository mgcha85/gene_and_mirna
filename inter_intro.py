import os
import pandas as pd
import sqlite3
import socket
from Database import Database
import sys
from Server import Server
from joblib import Parallel, delayed
import multiprocessing


class inter_intro:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def to_server(self):
        server = Server()
        server.connect(which='stokes')

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server.job_script(fname, time='00:30:00', which='stokes')

        server_root = os.path.join(server.server, 'source/gene_and_mirna')
        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def processInput(self, root, tname):
        print(tname)
        ref_fpath = os.path.join(self.root, 'software/stringtie-1.3.3b', 'genes.db')
        con_ref = sqlite3.connect(ref_fpath)

        fpath = os.path.join(root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.db')
        con_rsc = sqlite3.connect(fpath)

        out_path = fpath.replace('.db', '_out.db')
        out_con = sqlite3.connect(out_path)

        df_rsc = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con_rsc)

        for idx in df_rsc.index:
            chromosome = df_rsc.loc[idx, 'chromosome']
            start = df_rsc.loc[idx, 'start']
            end = df_rsc.loc[idx, 'end']
            strand = df_rsc.loc[idx, 'strand']

            tname_ref = 'genes_{}_{}'.format(chromosome, strand)
            df_ref = pd.read_sql_query("SELECT * FROM '{}' WHERE NOT start>{end} AND NOT end<{start} AND "
                                       "feature='transcript'".format(tname_ref, start=start, end=end), con_ref)
            if not df_ref.empty:
                df_rsc.loc[idx, 'type'] = 'intronic'
            else:
                df_rsc.loc[idx, 'type'] = 'intergenic'
        df_rsc.to_sql(tname, out_con, if_exists='replace', index=None)

    def get_reference(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.db')
        con_rsc = sqlite3.connect(fpath)

        tlist = Database.load_tableList(con_rsc)
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(self.processInput)(self.root, tname) for tname in tlist)


if __name__ == '__main__':
    ii = inter_intro()
    if ii.hostname == 'mingyu-Precision-Tower-7810':
        ii.to_server()
    else:
        ii.get_reference()
