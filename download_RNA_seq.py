import urllib
from bs4 import BeautifulSoup
import socket
import os
from fa_to_bed import fa2bed
import pandas as pd
import sqlite3


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

        server.job_script(fname, src_root=server_root, time='20:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
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
            for col in ["even col_28"]:
                items = sorted(soup.findAll("td", {"class": col}))
                for item in items:
                    contents.append(item.findAll("a")[0].attrs['href'])
                    download_url = item.findAll("a")[0].attrs['href']
                    print('download {}...'.format(download_url))
                    ulr_dir, fname = os.path.split(download_url)
                    urllib.request.urlretrieve(download_url, os.path.join(self.rna_dir, fname))

        with open('temp.csv', 'wt') as f:
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

    def to_bed(self):
        from time import time

        contents = self.get_file_pair()
        for fid, fpath in contents.items():
            start = time()
            sam_path = os.path.join(fid + '.sam')
            ret = self.f2b.comp_fa_to_sam(fpath, sam_path)
            if ret == 0:
                self.f2b.sam_to_bam(sam_path)
                self.f2b.bam_to_gtf(sam_path.replace('.bam', '.gtf'))
                # self.f2b.gff_to_gtf(sam_path.replace('.sam', '.gff'))
                # self.f2b.gtf_to_db(sam_path.replace('.sam', '.gtf'))

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
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            df = df[df['chromosome'].str.len() <= 5]
            df.to_sql(tname, con_out)


if __name__ == '__main__':
    drs = Download_RNA_seq()
    if drs.hostname == 'mingyu-Precision-Tower-7810':
        drs.to_server()
    else:
        drs.run()
