import pandas as pd
import os
import subprocess
import socket
from Server import Server
import sys
import shutil
import sqlite3
import numpy as np


class Convert:
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

        stdin, stdout, stderr = server.ssh.exec_command(
            "cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def command_exe(self, command):
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        return out

    def bam_to_gtf(self, bam_file):
        string_tie_root = os.path.join(self.root, 'software/stringtie-1.3.6')

        print('bam to gtf...')
        dirname, fname = os.path.split(bam_file)
        fname = os.path.splitext(fname)[0] + '.gtf'
        gtf_file = os.path.join(dirname, fname)

        command = '{str_root}/./stringtie -p 8 -G {str_root}/genes.gtf -o {gtf} -i {bam}' \
                  ''.format(str_root=string_tie_root, gtf=gtf_file, bam=bam_file)
        print(command)
        self.command_exe(command)

        print('done with bam to gtf')

    def gff_to_gtf(self, gff_file):
        string_tie_root = os.path.join(self.root, 'software/stringtie-1.3.6')

        print('gff to gtf...')
        dirname, fname = os.path.split(gff_file)
        fname = os.path.splitext(fname)[0] + '.gtf'
        gtf_file = os.path.join(dirname, fname)

        command = '{str_root}/./stringtie --merge -e -p 8 -G {str_root}/genes.gtf -o {gtf} -i {gff}' \
                  ''.format(str_root=string_tie_root, gtf=gtf_file, gff=gff_file)
        print(command)
        self.command_exe(command)
        print('done with bam to gtf')

    def gtf_to_db(self, gtf_file, con):
        gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        dirname, fname = os.path.split(gtf_file)

        df = pd.read_csv(gtf_file, names=gtf_columns, comment='#', sep='\t')
        df = df[df['feature'] == 'transcript'].reset_index(drop=True)
        df['chromosome'] = 'chr' + df['chromosome'].astype(str)
        attribute = df['attribute'].str.split('; ')

        pkm = []
        for idx in attribute.index:
            attr = attribute[idx]
            dict = {}
            for a in attr:
                key, value = a.split(' ')
                dict[key] = value.replace('"', '')

            row = [idx, None, None, dict['cov'], dict['FPKM'], dict['TPM'].replace(';', '')]
            if 'ref_gene_name' in dict:
                row[1] = dict['ref_gene_name']
                row[2] = dict['transcript_id']
            pkm.append(row)

        df_res = pd.DataFrame(data=pkm, columns=['index', 'gene_name', 'transcript_id', 'cov', 'FPKM', 'TPM'])
        df_res = df_res.set_index('index', drop=True)
        df = pd.concat([df, df_res], axis=1)
        df.drop(['attribute', 'frame'], axis=1).to_sql(os.path.splitext(fname)[0], con, index=None)

    def avg_rna_seq_by_tissues(self):
        df = pd.read_excel('RNA-seq_data_structure.xlsx')
        columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'gene_name', 'transcript_id']

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq.db')
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_tissue.db'))
        for src, df_src in df.groupby('tissue'):
            dfs = []
            index = []
            for fid in df_src['fid']:
                df_fid = pd.read_sql_query("SELECT * FROM '{}' WHERE FPKM>''".format(fid), con)
                df_fid.index = df_fid['chromosome'] + ':' + df_fid['start'].astype(str) + '-' + df_fid['end'].astype(
                    str)
                index.append(set(df_fid.index))
                dfs.append(df_fid)

            index = sorted(list(set.union(*index)))
            df_res = pd.DataFrame(index=index, columns=df_fid.columns)
            df_res[['cov', 'FPKM', 'TPM']] = 0
            for df_fid in dfs:
                df_res.loc[df_fid.index, ['cov', 'FPKM', 'TPM']] += df_fid[['cov', 'FPKM', 'TPM']].astype(float)
                df_res.loc[df_fid.index, columns] = df_fid[columns]

            df_res[['cov', 'FPKM', 'TPM']] /= len(dfs)
            df_res.to_sql(src, con_out, if_exists='replace', index=None)

    def stats_by_tissue(self):
        df = pd.read_excel('RNA-seq_data_structure.xlsx')
        contents = []
        df_grp = df.groupby('tissue')
        for tis, df_tis in df_grp:
            sources = ';'.join(df_tis['source'])
            fids = ';'.join(df_tis['fid'])
            contents.append([df_tis.shape[0], sources, fids])
        df_res = pd.DataFrame(contents, index=df_grp.groups, columns=['#data', 'source', 'fid'])
        df_res.to_excel('RNA-temp.xlsx')

    def arrange_files(self, bam_file):
        gtf_file = bam_file.replace('.bam', '.gtf')

        dirname, fname__ = os.path.split(bam_file)
        fname, ext = os.path.splitext(fname__)
        super_dir = '/'.join(dirname.split('/')[:-2])

        for from_path, type in zip([bam_file, gtf_file], ['bam', 'gtf']):
            to_path = os.path.join(super_dir, type, fname + '.' + type)
            shutil.move(from_path, to_path)

    def split_files(self, root_dir, batch_size, ext='.fastq'):
        flist = os.listdir(root_dir)
        flist = [os.path.join(root_dir, x) for x in flist if x.endswith(ext)]

        N = len(flist)
        M = (N + (batch_size - 1)) // batch_size
        for i, fpath in enumerate(flist):
            dirnum = (i // M) + 1
            _, fname = os.path.split(fpath)

            dirname = os.path.join(root_dir, str(dirnum))
            if not os.path.exists(dirname):
                os.mkdir(dirname)
            shutil.move(fpath, os.path.join(dirname, fname))

    def gtf_to_db_all(self, dirname):
        flist = [x for x in os.listdir(dirname + '/gtf') if x.endswith('.gtf')]
        con = sqlite3.connect(os.path.join(dirname, 'out', 'RNA_seq.db'))
        for fname in flist:
            self.gtf_to_db(os.path.join(dirname, 'gtf', fname), con)

    def run(self):
        # dirname = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/RNA-seq/bam/1'
        dirname = os.getcwd()

        flist = os.listdir(dirname)
        fpaths = [os.path.join(dirname, x) for x in sorted(flist) if x.endswith('.bam')]
        for fpath in fpaths:
            self.bam_to_gtf(fpath)
            self.arrange_files(fpath)


if __name__ == '__main__':
    con = Convert()
    if con.hostname == 'mingyu-Precision-Tower-7810':
        con.to_server()

    else:
        # con.stats_by_tissue()
        dirname = os.path.join(con.root, 'database', 'RNA-seq')
        con.gtf_to_db_all(dirname)

        # dirname = os.path.join(con.root, 'database', 'RNA-seq', 'bam')
        # flist = os.listdir(dirname)
        # for fname in flist:
        #     fpath = os.path.join(dirname, fname)
        #     con.bam_to_gtf(fpath)
        #     con.gtf_to_db(fpath.replace('.bam', '.gtf'))
