import pandas as pd
import numpy as np
import sqlite3
import os
import socket
import sys

from Database import Database
from Server import Server


class Correlation2:
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
        self.cells = []

    def to_server(self):
        server = Server()
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server.job_script(fname, time='03:00:00')

        server_root = os.path.join(server.server, 'source/gene_and_mirna')
        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def get_reference(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql_query("SELECT * FROM 'human_genes'", con)
        df_unique = df.drop_duplicates(subset=['gene_name'], keep=False)

        idx_duplicates = list(set(df.index) - set(df_unique.index))
        df_duplicates = df.loc[idx_duplicates]
        df_duplicates = df_duplicates[df_duplicates['source'] == 'ENSEMBL']
        df_duplicates = df_duplicates.drop_duplicates(subset=['gene_name'], keep=False)
        return pd.concat([df_unique, df_duplicates]).reset_index(drop=True)

    def pre_data(self):
        fpath = os.path.join(self.root, 'software/stringtie-1.3.6', 'genes.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        out_con = sqlite3.connect(fpath.replace('.db', '_2.db'))

        for tname in tlist:
            df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            attribute = df_rna['attribute'].str.split('; ')

            pkm = []
            for attr in attribute:
                dict = {}
                for a in attr[:-1]:
                    key, value = a.split(' ')
                    dict[key] = value.replace('"', '')

                if 'gene_name' not in dict:
                    pkm.append(np.nan)
                    continue
                pkm.append(dict['gene_name'])

            df_res = pd.DataFrame(data=pkm, columns=['gene_name'])
            df_rna = pd.concat([df_rna, df_res], axis=1)
            df_rna.drop(['attribute', 'frame'], axis=1).to_sql(tname, out_con, index=None, if_exists='append')

    def run(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con = sqlite3.connect(fpath)
        df_ref = pd.read_sql_query("SELECT * FROM 'human_genes_wo_duplicates'", con)
        SCOPE = 500

        fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con_rna = sqlite3.connect(fpath_rna)
        tlist_rna = Database.load_tableList(con_rna)

        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/tissues/out/', 'fantom_cage_by_tissue.db')
        con_fan = sqlite3.connect(fpath_fan)
        tlist_fan = Database.load_tableList(con_fan)

        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'correlation.db')
        out_con = sqlite3.connect(out_path)

        for tname in tlist_fan:
            contents = []
            tname_rna = [x for x in tlist_rna if tname in x]
            for rtname in tname_rna:
                print(rtname)
                for idx in df_ref.index:
                    chromosome = df_ref.loc[idx, 'chromosome']
                    strand = df_ref.loc[idx, 'strand']
                    if strand == '+':
                        tss = df_ref.loc[idx, 'start']
                    else:
                        tss = df_ref.loc[idx, 'end']

                    start = tss - SCOPE
                    end = tss + SCOPE
                    gene_name = df_ref.loc[idx, 'gene_name']

                    df_fan = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND NOT start>{end} AND "
                                               "NOT end<{start}".format(tname, chromosome, start=start, end=end), con_fan)
                    if df_fan.empty:
                        continue
                    counts = df_fan['score'].sum()

                    df_rna = pd.read_sql_query("SELECT * FROM '{}' WHERE gene_name='{}'".format(rtname, gene_name), con_rna)
                    if df_rna.empty:
                        continue

                    df_rna[['TPM', 'FPKM']] = df_rna[['TPM', 'FPKM']].astype(float)
                    tpm = df_rna['TPM'].mean()
                    fpkm = df_rna['FPKM'].mean()
                    loci = ';'.join([chromosome, tss, strand])
                    contents.append([gene_name, loci, tpm, fpkm, counts, rtname])

                df_rep = pd.DataFrame(contents, columns=['gene_name', 'loci', 'tpm', 'fpkm', 'counts', 'tissue'])
                df_rep.to_sql(rtname, out_con, if_exists='replace', index=None)


if __name__ == '__main__':
    cor = Correlation2()
    if cor.hostname == 'mingyu-Precision-Tower-7810':
        # cor.to_server()
        cor.run()
    else:
        cor.run()
