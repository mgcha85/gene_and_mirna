import urllib
from bs4 import BeautifulSoup
import socket
import os
from fa_to_bed import fa2bed
import pandas as pd
import sqlite3


class Download_RNA_seq:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.cage_dir = os.path.join(self.root, 'database/Fantom/v5/tissues')
        self.url = 'http://fantom.gsc.riken.jp/5/sstar/'
        self.f2b = fa2bed()

    def get_script(self, url):
        try:
            with urllib.request.urlopen(url) as response:
                html = response.read().decode('utf-8')
                return html
        except Exception as e:
            raise e

    def get_cage_peak_id(self):
        fpath = os.path.join(self.root, 'Papers/complement', 'supp_gkv608_nar-00656-h-2015-File009.xls')
        df = pd.read_excel(fpath, header=1)
        return df.dropna(axis=0, how='any')['CAGE FF_ID']

    def to_db(self, fpath):
        self.f2b.bam_to_gtf(fpath)
        self.f2b.gtf_to_db(fpath.replace('.bam', '.gtf'))

    def pc_run(self, cid, i, N):
        print('{:.2f}%'.format(100 * (i + 1) / N))

        url = self.url + 'FF:' + cid
        page = self.get_script(url)
        soup = BeautifulSoup(page, "lxml")
        items = soup.findAll("a", {"class": 'external text'})

        content = None
        for item in items:
            download_url = item.attrs['href']
            if '{}.hg19.nobarcode.bam'.format(cid) in download_url:
                ulr_dir, fname = os.path.split(download_url)
                tissue = fname.split('%')[0]
                local_dir = os.path.join(self.cage_dir, tissue)
                if not os.path.exists(local_dir):
                    os.mkdir(local_dir)

                local_path = os.path.join(local_dir, fname)
                content = [download_url, local_path]
                urllib.request.urlretrieve(download_url, local_path)
                self.to_db(local_path)
        return content

    def run(self):
        cage_id = self.get_cage_peak_id()[2:]
        N = len(cage_id)

        contents = []
        for i, cid in enumerate(cage_id):
            contents.append(self.pc_run(cid, i, N))

        df = pd.DataFrame(data=contents, columns=['download_url', 'loca_path'])
        df.to_csv('download_cage_peak_id.csv', index=None)

    def check_not_download(self):
        with open('temp.csv') as f:
            down_list__ = f.read().split('\n')

        down_list = {}
        for dl in down_list__:
            dirname, fname = os.path.split(dl)
            down_list[fname] = dirname

        dirname = os.path.join(self.root, 'database', 'RNA-seq')
        flist = os.listdir(dirname)
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

    def to_bed(self):
        from time import time

        contents = self.get_file_pair()
        for fid, fpath in contents.items():
            start = time()
            sam_path = os.path.join(fid + '.sam')
            ret = self.f2b.comp_fa_to_sam(fpath, sam_path)
            if ret == 0:
                self.f2b.sam_to_bam(sam_path)
                self.f2b.bam_to_gtf(sam_path.replace('.sam', '.bam'))

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
    drs.run()