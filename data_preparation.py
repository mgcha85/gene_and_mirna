import os
import pandas as pd
import numpy as np
import sqlite3
import socket
from Database import Database
import sys
from Server import Server


class data_preparation:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinfomatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        self.bed_columns = ['chromosome', 'start', 'end', 'location', 'score', 'strand']

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

    # 1. data download from FANTOM5
    def download_cell_line_fantom(self):
        import urllib.request

        df = pd.read_csv('fantom_cell_line_list.txt')
        for idx in df.index:
            link = df.loc[idx, 'link']
            url, fname = os.path.split(link)
            dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
            fpath = os.path.join(dirname, fname)
            urllib.request.urlretrieve(link, fpath)
            # df.loc[idx, 'link'] = link.replace('%', '%25')
        # df.to_csv('fantom_cell_line_list.txt', index=None)

    # 2. check if it is hg19 or not
    def check_chain(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        fpath = os.path.join(dirname, 'list/00_human.cell_line.hCAGE.hg19.assay_sdrf2.xlsx')
        df_list = pd.read_excel(fpath)
        df_sub = df_list[~df_list['File Name.1'].str.contains('hg19')]
        print(df_sub)

    def merge_processInput(self, row, columns):
        fpath, cline, ename = row
        return pd.read_csv(fpath, sep='\t', compression='gzip', names=columns), cline, ename

    # 3. integrate cell line data
    def merge_cline_db(self):
        from joblib import Parallel, delayed
        num_cores = 6

        def as_batch(data, i, batch_size=100):
            N = len(data)
            sidx = i * batch_size
            eidx = sidx + batch_size
            if eidx > N:
                eidx = N
            return data[sidx: eidx]

        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        fpath = os.path.join(dirname, 'list/00_human.cell_line.hCAGE.hg19.assay_sdrf2.xlsx')
        df_list = pd.read_excel(fpath)
        df_grp = df_list.groupby('cell line')
        columns = ['chromosome', 'start', 'end', 'location', 'score', 'strand']

        con = sqlite3.connect(os.path.join(dirname, 'human_cell_line_hCAGE.db'))
        fpaths = []
        fails = []
        multi = []
        for cline, df_sub in df_grp:
            if df_sub.shape[0] > 1:
                print(cline, df_sub.shape[0])
                multi.append([cline, df_sub.shape[0]])
            idx = df_sub.index[0]
            fname = df_sub.loc[idx, 'File Name.1']
            ename = df_sub.loc[idx, 'Extract Name']
            fname = fname.replace('.nobarcode.bam', '.ctss.bed.gz')

            fpath = os.path.join(dirname, fname)
            if not os.path.exists(fpath):
                fails.append(fname)
                continue
            fpaths.append([fpath, cline, ename])

        with open('fails.txt', 'wt') as f:
            f.write('\n'.join(fails))
        df_multi = pd.DataFrame(multi, columns=['cell line', 'N'])
        df_multi.to_csv('multiple_cell_lines.csv', index=None)

        N = int((len(fpaths) + num_cores) // num_cores)
        for i in range(N):
            dfs = Parallel(n_jobs=num_cores)(delayed(self.merge_processInput)(row, columns) for row in as_batch(fpaths, i, num_cores))
            for row in dfs:
                df, cline, ename = row
                df[['chromosome', 'start', 'end', 'strand', 'score']].to_sql(cline, con, if_exists='replace', index=None)


if __name__ == '__main__':
    dp = data_preparation()
    if dp.hostname == 'mingyu-Precision-Tower-7810':
        # dp.db_to_bed()
        dp.download_cell_line_fantom()
        # dp.bed_to_db()
    else:
        dp.bed_to_db()

