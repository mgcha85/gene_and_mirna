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
            self.root = '/home/mingyu/Bioinformatics'
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

    def download_tissue_fantom(self):
        import urllib.request

        url__ = 'http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.tissue.hCAGE/{}'
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        fpath = os.path.join(dirname, '00_human.tissue.hCAGE.hg19.assay_sdrf.csv')
        df = pd.read_csv(fpath, sep='\t')
        df['Comment [sample_name]'] = df['Comment [sample_name]'].str.split(', ').str[0]
        df['Comment [sample_name]'] = df['Comment [sample_name]'].str.split(' - ').str[0]

        for grp, df_grp in df.groupby('Comment [sample_name]'):
            subdir = os.path.join(dirname, grp)
            if not os.path.exists(subdir):
                os.mkdir(subdir)
            for idx in df_grp.index:
                fname = df_grp.loc[idx, 'File Name.1'].replace('.nobarcode.bam', '.ctss.bed.gz').replace('%', '%25')
                print(grp, fname)
                url = url__.format(fname)
                urllib.request.urlretrieve(url, os.path.join(subdir, fname))

    def correlation_replicates(self, df, cols):
        from itertools import combinations

        report = pd.Series()
        for comb in combinations(cols, 2):
            coeff = np.corrcoef(df[comb[0]].fillna(0), df[comb[1]].fillna(0))
            report['&'.join(comb)] = coeff[0, 1]
        return report

    def avg_fantom_by_tissue(self):
        columns = ['chromosome', 'start', 'end', 'loc', 'score', 'strand']
        dirname_ = os.path.join(self.root, 'database/Fantom/v5/tissues')

        fpath = os.path.join(dirname_, 'FANTOM_tissue.db')
        con_out = sqlite3.connect(fpath)

        # writer = pd.ExcelWriter(os.path.join(self.root, 'database/Fantom/v5/tissues', 'duplicates.xlsx'), engine='xlsxwriter')
        flist = [os.path.join(dirname_, x) for x in os.listdir(dirname_)]
        dirlist = sorted([x for x in flist if os.path.isdir(x)])

        for dirname in dirlist[72:]:
            tissue = dirname.split('/')[-1]
            dfs = []
            index = []

            flist = os.listdir(dirname)
            flist = [f for f in flist if f.endswith('.gz')]
            print(tissue, len(flist))

            for fname in flist:
                fpath = os.path.join(dirname, fname)
                df_fid = pd.read_csv(fpath, compression='gzip', sep='\t', names=columns)
                df_fid.index = df_fid['chromosome'] + ':' + df_fid['start'].astype(str) + '-' + df_fid['end'].astype(str)
                index.append(set(df_fid.index))
                dfs.append(df_fid)

            index = sorted(list(set.union(*index)))
            df_res = pd.DataFrame(index=index, columns=df_fid.columns)
            cols = []
            for i, df_fid in enumerate(dfs):
                cols.append('score (Rep {})'.format(i + 1))
                df_res.loc[df_fid.index, cols[-1]] = df_fid['score'].astype(float)
                df_res.loc[df_fid.index, cols[-1]] /= df_res.loc[df_fid.index, cols[-1]].sum() / 1e6
                if i == 0:
                    df_res.loc[df_fid.index, columns] = df_fid[columns]

            # report = self.correlation_replicates(df_res, cols)
            # report.to_excel(writer, sheet_name=tissue)
            df_res['score'] = df_res[cols].mean(axis=1)
            print(df_res[cols].sum(axis=0))
            df_res.drop(['loc'], axis=1).to_sql(tissue, con_out, if_exists='replace', index=None)
        # writer.save()
        # writer.close()


if __name__ == '__main__':
    dp = data_preparation()
    if dp.hostname == 'mingyu-Precision-Tower-7810':
        # dp.db_to_bed()
        dp.to_server()
        # dp.bed_to_db()
    else:
        # dp.avg_fantom_by_tissue()

        from Correlation_gpu import Correlation
        from Regression import Regression

        cor = Correlation()
        rg = Regression()

        cor.sum_fan(100, ref='gene')
        cor.sum_fan(100, ref='mir')
        rg.regression(100, 'nz')
