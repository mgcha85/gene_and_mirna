import os
import sqlite3
import pandas as pd
import numpy as np
import socket
from Database import Database
from Server import Server
import sys


class Fantom_RNA:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.cells = []
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

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

    def get_tissues(self):
        df = pd.read_csv("E-MTAB-1733.csv")
        tissues = np.array([x.split('_')[0] for x in df['src']])
        tissues = np.unique(tissues)
        with open('tissues.txt', 'wt') as f:
            f.write('\n'.join(tissues))

    def extract_tissues_from_fantom(self):
        with open('tissues.txt', 'rt') as f:
            tissues = f.read().split('\n')

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df = pd.read_csv(fpath)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.tissues.db')
        con = sqlite3.connect(fpath)

        for tissue in tissues:
            print(tissue)
            columns = ['chromosome', 'start', 'end', 'strand', 'length']
            for col in df.columns[5:]:
                if tissue in col:
                    columns.append(col)
            df[columns].to_sql(tissue, con, if_exists='replace', index=None)

    def tpm(self, df_gene, df_rna):
        length = abs(df_gene['start'] - df_gene['end'])
        N = df_rna.shape[0]
        rpk = N / length
        return rpk / 1e6

    def rpkm(self, df_gene, df_rna):
        length = abs(df_gene['start'] - df_gene['end'])
        N = df_rna.shape[0]
        rpm = N / 1e6
        return rpm / length

    def merge_db(self):
        df = pd.read_csv("E-MTAB-1733.csv")
        new_col = [''] * df.shape[0]
        for i, url in enumerate(df.link):
            addr, fname = os.path.split(url)
            fid = os.path.splitext(fname)[0].split('_')[0]
            new_col[i] = fid
        df['fid'] = new_col
        df = df.drop_duplicates(subset=['fid'])
        df = df.set_index('fid', drop=True)

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_full.db')
        con = sqlite3.connect(fpath)

        dirname = os.path.join(self.root, 'database/RNA-seq/gtf')
        flist = os.listdir(dirname)
        N = len(flist)
        for i, fid in enumerate(flist):
            print('{} / {}'.format(i + 1, N))
            fid, ext = os.path.splitext(fid)
            tissue = df.loc[fid, 'src']

            df_sub = pd.read_csv(os.path.join(dirname, fid + ext), names=self.gtf_columns, sep='\t', comment='#')
            # df_sub = df_sub[df_sub['feature'] == 'transcript'].reset_index(drop=True)
            df_sub['chromosome'] = 'chr' + df_sub['chromosome'].astype(str)
            attribute = df_sub['attribute'].str.split('; ')

            pkm = []
            for attr in attribute:
                dict = {}
                for a in attr:
                    key, value = a.split(' ')
                    dict[key] = value.replace('"', '')

                if 'ref_gene_name' not in dict:
                    pkm.append([np.nan] * 2)
                    continue
                pkm.append([dict['ref_gene_name'], dict['cov']])

            df_res = pd.DataFrame(data=pkm, columns=['gene_name', 'cov'])
            df_sub = pd.concat([df_sub, df_res], axis=1)
            df_sub.drop(['attribute', 'frame'], axis=1).to_sql(tissue, con, index=None, if_exists='append')

    def processInput(self, df, tissue, idx):
        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq.db')
        con_rna = sqlite3.connect(fpath)

        if idx % 1000 == 0 or idx + 1 == df.shape[0]:
            print('{:0.2f}%'.format(100 * (idx + 1) / df.shape[0]))

        chromosome = df.loc[idx, 'chromosome']
        start = df.loc[idx, 'start']
        end = df.loc[idx, 'end']
        strand = df.loc[idx, 'strand']

        df_rna = pd.read_sql("SELECT * FROM '{}_4a' WHERE chromosome='{}' AND strand='{}' AND NOT start>{end} "
                                   "AND NOT stop<{start}".format(tissue, chromosome, strand, start=start, end=end),
                                   con_rna)
        return self.rpkm(df.loc[idx], df_rna)

    def check_empty(self):
        df = pd.read_csv('E-MTAB-1733.csv')
        fid = []
        for idx in df.index:
            url = df.loc[idx, 'link']
            _, fname = os.path.split(url)
            fname, _ = os.path.splitext(fname)
            fname = fname.split('_')[0]
            fid.append(fname)
        fid = np.unique(np.array(fid))

        dirname = os.path.join(self.root, 'database/RNA-seq')
        flist = os.listdir(dirname)
        flist = [os.path.splitext(x)[0] for x in flist if x.endswith('.db')]

        remains = list(set(fid) - set(flist))
        print(remains)

    def run(self):
        from Database import Database
        from joblib import Parallel, delayed
        import multiprocessing
        # num_cores = 6
        num_cores = multiprocessing.cpu_count() // 2
        print('num_cores: {}'.format(num_cores))

        fpath_out = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.out.db')
        con_out = sqlite3.connect(fpath_out)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.tissues.db')
        con = sqlite3.connect(fpath)
        tissues = Database.load_tableList(con)

        for tissue in tissues:
            print(tissue)
            df = pd.read_sql("SELECT * FROM '{}'".format(tissue), con)
            if df.shape[1] <= 5:
                continue

            df['RPKM'] = Parallel(n_jobs=num_cores)(delayed(self.processInput)(df, tissue, idx) for idx in df.index)
            df.to_sql(tissue, con_out, index=None, if_exists='replace')


if __name__ == '__main__':
    fr = Fantom_RNA()
    if fr.hostname == 'mingyu-Precision-Tower-7810':
        # fr.merge_db()
        fr.to_server()
    else:
        fr.merge_db()
