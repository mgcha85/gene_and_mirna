import pandas as pd
import numpy as np
import sqlite3
import os
import socket
from copy import deepcopy
from Database import Database


class data_eval:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

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

    def run(self):
        dirname = os.path.join(self.root, 'amlan')
        flist = [x for x in os.listdir(dirname) if x.endswith('.txt')]

        # load all data
        dfs = []
        for fname in flist:
            fpath = os.path.join(dirname, fname)
            tissue = '_'.join(os.path.splitext(fname)[0].split('_')[4:])
            print('load {}...'.format(tissue))
            df = pd.read_csv(fpath, sep='\t', index_col=2)
            # df.index.name = 'transcript_id'
            df = df.rename(columns={'HPA RNA-Seq': 'HPA RNA-Seq:{}'.format(tissue),
                                    'CAGE': 'CAGE:{}'.format(tissue)})
            dfs.append(df)

        df_res = pd.concat(dfs, axis=1)
        # fan_cols = [col.split(':')[1] for col in df_res.columns if 'CAGE' in col]
        fan_cols = [col for col in df_res.columns if 'CAGE' in col]
        rna_cols = [col for col in df_res.columns if 'HPA RNA-Seq' in col]

        con = sqlite3.connect(os.path.join(dirname, 'gene_exp_wit_cage_tissue.db'))
        df_fan = df_res[fan_cols]
        df_rna = df_res[rna_cols]

        df_fan.columns = [col.split(':')[1] for col in df_rna.columns]
        df_rna.columns = [col.split(':')[1] for col in df_rna.columns]

        df_fan.to_sql('CAGE', con, if_exists='replace')
        df_rna.to_sql('RNA_seq', con, if_exists='replace')

    def add_location(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v19.annotation.gtf.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT transcript_id, chromosome, start, end, strand FROM 'gencode.v19.annotation.gtf' GROUP BY transcript_id", con, index_col='transcript_id')

        dirname = os.path.join(self.root, 'amlan')
        con = sqlite3.connect(os.path.join(dirname, 'gene_exp_wit_cage_tissue.db'))
        for lab in ['CAGE', 'RNA_seq']:
            df_tis = pd.read_sql("SELECT * FROM '{}'".format(lab), con, index_col='Transcript ID')
            df_tis = pd.concat([df_tis, df.loc[df_tis.index]], axis=1)
            df_tis.loc[df_tis.index].to_sql('{}_loc'.format(lab), con, if_exists='replace')

    def check_non_one(self):
        dirname = os.path.join(self.root, 'amlan/corr')
        fpath = os.path.join(dirname, 'correlated_transcripts.txt')
        df = pd.read_csv(fpath, sep='\t', names=['transcript_id', 'corr']).set_index('transcript_id')
        df = df[df['corr'] < 1]
        print(df.shape[0])

    def scan(self):
        dirname = os.path.join(self.root, 'amlan')
        fpath = os.path.join(dirname, 'corr', 'correlated_transcripts.txt')
        df = pd.read_csv(fpath, sep='\t', names=['transcript_id', 'corr']).set_index('transcript_id')
        df = df[df['corr'] == 1]

        con = sqlite3.connect(os.path.join(dirname, 'gene_exp_wit_cage_tissue.db'))
        df_rna = pd.read_sql("SELECT * FROM 'RNA_seq'", con, index_col='Transcript ID')
        df_cage = pd.read_sql("SELECT * FROM 'CAGE'", con, index_col='Transcript ID')

        df_rna = df_rna.loc[df.index]
        df_cage = df_cage.loc[df.index]
        print('see')

    def comparison(self):
        df__ = []
        tid = []

        for label in ['_raw', '_processed', '_rand']:
            fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100{}.db'.format(label))
            con = sqlite3.connect(fpath)
            dfs = []
            for tname in Database.load_tableList(con):
                dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='transcript_id'))
            df = pd.concat(dfs)
            tid.append(set(df.index))
            mean = df['corr'].mean()
            median = df['corr'].median()
            std = df['corr'].std()
            min = df['corr'].min()
            max = df['corr'].max()
            q7 = df['corr'].quantile(0.7)
            q8 = df['corr'].quantile(0.8)
            q9 = df['corr'].quantile(0.9)
            print('[{}] mean: {:0.2f}, median: {:0.2f}, std: {:0.2f}, min: {:0.2f}, max: {:0.2f}, q7: {:0.2f}, '
                  'q8: {:0.2f}, q9:{:0.2f}'.format(label, mean, median, std, min, max, q7, q8, q9))
            df__.append(df['corr'])

        tid = set.intersection(*tid)
        print(len(tid))
        df_res = pd.concat(df__, axis=1)
        df_res.columns = ['raw', 'processed', 'rand']
        df_res = df_res.loc[tid]
        corr = df_res.corr(method='spearman')
        print(corr)

    def random_gen(self, type='rna-seq'):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_100_raw.db')
        con = sqlite3.connect(fpath)
        tlist = [x for x in Database.load_tableList(con) if type in x]

        dfs = []
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='transcript_id')
            dfs.append(df)
        df = pd.concat(dfs)
        values = df.loc[:, 'appendix':]
        values = values.replace(-1, 0)
        max, min = values.values.max(), values.values.min()
        diff = max - min
        rand = np.random.random(size=(df.shape[0], 22))
        rand = rand * diff + min

        out_con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_100_rand.db'))
        pd.DataFrame(rand, columns=values.columns, index=values.index).to_sql(type, out_con, if_exists='replace')

    def random_corr(self):
        con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_100_rand.db'))
        df_rna = pd.read_sql("SELECT * FROM 'rna-seq'", con, index_col='transcript_id')
        df_fan = pd.read_sql("SELECT * FROM 'fantom'", con, index_col='transcript_id')

        from corr_gpu import Spearman
        con_out = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100_rand.db'))
        corr = Spearman(self.root)
        df_res = pd.Series(data=corr.run(df_fan, df_rna, prod=False), index=df_fan.index)
        df_res.to_sql('corr', con_out, if_exists='replace')


if __name__ == '__main__':
    de = data_eval()
    if de.hostname == 'mingyu-Precision-Tower-781':
        de.to_server()
    else:
        # de.run()
        # de.random_gen()
        # de.random_corr()
        de.comparison()
        # de.check_non_one()
        # de.scan()
