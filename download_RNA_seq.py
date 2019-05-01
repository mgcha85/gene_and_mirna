import urllib
from bs4 import BeautifulSoup
import socket
import os
from fa_to_bed import fa2bed
import pandas as pd


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
        self.rna_dir = os.path.join(self.root, 'database/RNA-seq/9')
        self.url = 'https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1733/samples/'
        self.f2b = fa2bed()

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
        page_size = 25
        iter = int((N + page_size - 1) // page_size)

        contents = []
        for i in range(1, iter):
            url = self.url + '?s_page={}&s_pagesize=25'.format(i)
            page = self.get_script(url)

            soup = BeautifulSoup(page, "lxml")
            for col in ["odd col_28", "even col_28"]:
                items = soup.findAll("td", {"class": col})
                for item in items:
                    contents.append(item.findAll("a")[0].attrs['href'])
                    download_url = item.findAll("a")[0].attrs['href']
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
        fileList = os.listdir(self.rna_dir)
        contents = {}
        for fname in fileList:
            if ext not in fname:
                continue
            fpath = os.path.join(self.rna_dir, fname)
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
            sam_path = os.path.join(self.rna_dir, fid + '.sam')
            ret = self.f2b.comp_fa_to_sam(fpath, sam_path)
            if ret == 0:
                self.f2b.sam_to_bam(sam_path)
                self.f2b.bam_to_gtf(sam_path.replace('.sam', '.bam'))

            elapsed = time() - start
            print('[{}] {:.0f}:{}'.format(fid, elapsed // 60, elapsed % 60))


if __name__ == '__main__':
    drs = Download_RNA_seq()
    drs.to_bed()
