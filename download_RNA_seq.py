import urllib
from bs4 import BeautifulSoup
import socket
import os
from fa_to_bed import fa2bed
import pandas as pd
import sqlite3
from collections import OrderedDict
from Convert import Convert
from Util import Util
from Database import Database


class Download_RNA_seq:
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
        self.rna_dir = os.path.join(self.root, 'database/RNA-seq/fastq')
        self.url = 'https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1733/samples/'
        self.f2b = fa2bed()

    def to_server(self, tag):
        from Server import Server
        import sys

        which = 'stokes'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname__ = os.path.split(local_path)
        fname = fname__.replace('.py', '_{}.py'.format(tag))

        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)
        server_path = server_path.replace(fname__, fname)

        server.job_script(fname, src_root=server_root, time='24:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command(
            "cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def sort_to_done(self):
        import shutil

        root = os.path.join(self.root, 'database/RNA-seq')
        # flist = os.listdir(root)
        # flist = [os.path.splitext(x)[0] for x in flist if x.endswith('.db')]

        with open('flist.txt', 'rt') as f:
            flist = f.read().split('\n')

        sub_flists = []
        for i in range(1, 9):
            sub_dirname = os.path.join(root, str(i))
            sub_flist = os.listdir(sub_dirname)
            sub_flist = [os.path.join(sub_dirname, x) for x in sub_flist]
            sub_flists.extend(sub_flist)

        contents = []
        for db_fname in flist:
            for sfname in sub_flists:
                if db_fname in sfname:
                    contents.append(sfname)

        contents = sorted(contents)
        for con in contents:
            dirname, fname = os.path.split(con)
            print('remove: ', con)
            os.remove(con)
            # out_path = os.path.join(root, 'done', fname)
            # print(con, ' --> ', out_path)
            # shutil.move(con, out_path)

    def get_script(self, url):
        try:
            with urllib.request.urlopen(url) as response:
                html = response.read().decode('utf-8')
                return html
        except Exception as e:
            raise e

    def get_header(self, page):
        soup = BeautifulSoup(page, "lxml")
        script = soup.findAll("div", {"class": "ae-stats"})
        for scr in script:
            contents = scr.findAll("span")[1]
            return int(contents.text)

    def run(self):
        page = self.get_script(self.url)
        N = self.get_header(page)
        page_size = 500
        iter = int((N + page_size - 1) // page_size)

        contents = []
        for i in range(1, iter + 1):
            url = self.url + '?s_page={}&s_pagesize={}'.format(i, page_size)
            page = self.get_script(url)

            soup = BeautifulSoup(page, "lxml")
            for col in ["even col_28", "odd col_28"]:
                items = soup.findAll("td", {"class": col})
                for item in items:
                    download_url = item.findAll("a")[0].attrs['href']
                    contents.append(download_url)
                    print('download {}...'.format(download_url))
                    # ulr_dir, fname = os.path.split(download_url)
                    # urllib.request.urlretrieve(download_url, os.path.join(self.rna_dir, fname))

        with open('RNA-seq_list.csv', 'wt') as f:
            f.write('\n'.join(sorted(contents)))

    def check_not_download(self):
        with open('temp.csv') as f:
            down_list__ = f.read().split('\n')
        
        down_list = {}
        for dl in down_list__:
            dirname, fname = os.path.split(dl)
            down_list[fname] = dirname

        dirname = os.path.join(self.root, 'database', 'RNA-seq/fastq')
        flist = []
        for path, subdirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.gz'):
                    flist.append(name)
        remains = list(set(down_list.keys()) - set(flist))

        N = len(remains)
        if N == 0:
            print('download completed!!')
            return

        for i, fname in enumerate(remains):
            print('{} / {} [{}]'.format(i + 1, N, fname))
            dirname = down_list[fname]
            download_url = os.path.join(dirname, fname)
            ulr_dir, fname = os.path.split(download_url)
            urllib.request.urlretrieve(download_url, os.path.join(self.rna_dir, fname))

    def get_file_pair(self, ext='.gz'):
        dirname = os.getcwd()
        fileList = os.listdir(dirname)
        print(fileList)
        contents = {}
        for fname in fileList:
            if ext not in fname:
                continue
            fpath = os.path.join(dirname, fname)
            fname, ext = os.path.splitext(fname)
            fid = fname.split('_')[0]
            if fid in contents:
                contents[fid].append(fpath)
            else:
                contents[fid] = [fpath]
        return contents

    def get_pair(self, contents):
        pair_list = OrderedDict()
        for url in contents:
            key = url.split('/')[-2]
            if key not in pair_list:
                pair_list[key] = [url]
            else:
                pair_list[key].append(url)
        return pair_list

    def bam_to_bed(self, dirname):
        from time import time

        flist = []
        for path, subdirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.bam'):
                    flist.append(name)

        for fname in flist:
            fid, ext = os.path.splitext(fname)
            bam_path = os.path.join(fid + '.bam')
            start = time()

            self.f2b.bam_to_gff(bam_path)
            self.f2b.gff_to_gtf(bam_path.replace('.bam', '.gff'))

            elapsed = time() - start
            print('[{}] {:.0f}:{:.0f}'.format(fid, elapsed // 60, elapsed % 60))

    def to_bed(self, i, batchsize=10, download=False):
        from time import time

        if download:
            with open('RNA-seq_list.csv', 'rt') as f:
                rna_seq = f.read().split('\n')
                pair_list = self.get_pair(rna_seq)

            fids = list(pair_list.keys())[i*batchsize:(i+1)*batchsize]
            for fid in fids:
                fpaths = []
                for url in pair_list[fid]:
                    dirname, fname = os.path.split(url)
                    fpath = os.path.join(self.rna_dir, fname)
                    fpaths.append(fpath)
                    urllib.request.urlretrieve(url, fpath)

                sam_path = os.path.join(self.rna_dir, fid + '.sam')
                self.f2b.comp_fa_to_sam(fpaths, sam_path)
                self.f2b.sam_to_bam(sam_path)
                self.f2b.bam_to_gtf(sam_path.replace('.sam', '.bam'))

        else:
            pair_list = self.get_file_pair()
            for fid, fpath in pair_list.items():
                start = time()
                sam_path = os.path.join(fid + '.sam')
                ret = self.f2b.comp_fa_to_sam(fpath, sam_path)
                if ret == 0:
                    self.f2b.sam_to_bam(sam_path)
                    self.f2b.bam_to_gtf(sam_path.replace('.bam', '.gtf'))
                    self.f2b.gtf_to_db(self.root)

            elapsed = time() - start
            print('[{}] {:.0f}:{}'.format(fid, elapsed // 60, elapsed % 60))

    def merge_rna_seq(self):
        dirname = os.path.join(self.root, 'database', 'RNA-seq')
        flist = os.listdir(dirname)
        flist = [x for x in flist if x.endswith('.db')]

        con_out = sqlite3.connect(os.path.join(dirname, 'out', 'RNA_seq.db'))
        for fname in flist:
            con = sqlite3.connect(os.path.join(dirname, fname))
            tname = os.path.splitext(fname)[0]
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            df = df[df['chromosome'].str.len() <= 5]
            df.to_sql(tname, con_out)

    def avg_by_rep(self):
        con = sqlite3.connect(os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq.db'))
        con_out = sqlite3.connect(os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_src.db'))
        df = pd.read_excel("RNA-seq_data_structure.xlsx")
        for src, df_sub in df.groupby('source'):
            print(src)
            dfs, tids, cols = [], [], []
            for i, fid in enumerate(df_sub['fid']):
                df_rna = pd.read_sql("SELECT chromosome, start, end, strand, ref_gene_id, ref_gene_name, "
                                     "reference_id, FPKM as FPKM_{fid} FROM '{fid}'".format(fid=fid), con,
                                     index_col='reference_id')
                dfs.append(df_rna)
                cols.append('FPKM_{}'.format(fid))
                tids.append(set(df_rna.index))

            columns__ = ['chromosome', 'start', 'end', 'strand', 'ref_gene_id', 'ref_gene_name']
            df_merge = pd.DataFrame(index=list(set.union(*tids)), columns=columns__ + cols)
            for i, df_sub in enumerate(dfs):
                columns = columns__ + [cols[i]]
                df_merge.loc[df_sub.index, columns] = df_sub.loc[:, columns]

            df_merge[cols] = df_merge[cols].astype(float)
            df_merge['FPKM'] = df_merge[cols].mean(axis=1)
            df_merge.index.name = 'reference_id'
            df_merge.to_sql(src, con_out, if_exists='replace')

    def avg_by_tissue(self):
        con = sqlite3.connect(os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_src.db'))
        con_out = sqlite3.connect(os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db'))
        tlist = Database.load_tableList(con)

        table_dict = {}
        for tname in tlist:
            tissue = tname.split('_')[0]
            if tissue not in table_dict:
                table_dict[tissue] = [tname]
            else:
                table_dict[tissue].append(tname)

        for tissue, tnames in table_dict.items():
            print(tissue)
            dfs, tids, cols = [], [], []
            for tname in tnames:
                df = pd.read_sql("SELECT chromosome, start, end, strand, ref_gene_id, ref_gene_name, "
                                         "reference_id, FPKM as FPKM_{tname} FROM '{tname}'"
                                 "".format(tname=tname), con, index_col='reference_id')
                dfs.append(df)
                cols.append('FPKM_{}'.format(tname))
                tids.append(set(df.index))

            columns__ = ['chromosome', 'start', 'end', 'strand', 'ref_gene_id', 'ref_gene_name']
            df_merge = pd.DataFrame(index=list(set.union(*tids)), columns=columns__ + cols)
            for i, df_sub in enumerate(dfs):
                columns = columns__ + [cols[i]]
                df_merge.loc[df_sub.index, columns] = df_sub.loc[:, columns]

            df_merge[cols] = df_merge[cols].astype(float)
            df_merge[cols] = df_merge[cols].fillna(0)
            df_merge['FPKM'] = df_merge[cols].mean(axis=1)
            df_merge.index.name = 'reference_id'
            df_merge.to_sql(tissue, con_out, if_exists='replace')

    def split(self):
        con = sqlite3.connect(os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db'))
        con_out = sqlite3.connect(os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_spt.db'))

        for tname in Database.load_tableList(con):
            print(tname)
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            for str, df_str in df.groupby('strand'):
                for chr, df_chr in df_str.groupby('chromosome'):
                    if len(chr) <= 5:
                        df_chr.drop(['chromosome', 'strand'], axis=1).to_sql('_'.join([tname, chr, str]), con_out, index=None, if_exists='replace')


if __name__ == '__main__':
    drs = Download_RNA_seq()
    batch_size = 10

    if drs.hostname == 'mingyu-Precision-Tower-781':
        drs.to_server("")
        # with open('RNA-seq_list.csv', 'rt') as f:
        #     rna_seq = f.read().split('\n')
        #     pair_list = drs.get_pair(rna_seq)
        #
        #     n = len(pair_list)
        #     m = int((n + batch_size - 1) / batch_size)
        #     for i in range(m):
        #         drs.to_server(str(i))
    else:
        # drs.avg_by_rep()
        drs.avg_by_tissue()
        drs.split()

        # drs.f2b.gtf_to_db(os.path.join(drs.root, 'database/RNA-seq'))

        # conv = Convert()
        # conv.avg_rna_seq_by_tissues()
        #
        # util = Util()
        # fpath = os.path.join(drs.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        # con = sqlite3.connect(fpath)
        # for tname in Database.load_tableList(con):
        #     util.split(fpath, tname)

        # import sys
        #
        # local_path = sys.argv[0]
        # dirname, fname = os.path.split(local_path)
        # fname_, ext = os.path.splitext(fname)
        # tag = fname_.split('_')[-1]
        # drs.to_bed(int(tag), batch_size, download=True)
        # drs.to_bed(0, batch_size, download=True)
