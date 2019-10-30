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
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
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

    def get_attr(self, attr):
        dfs = []
        for att in attr:
            at = [a.split(' ') for a in att]
            df = pd.DataFrame(at).set_index(0)
            df.iloc[:, -1] = df.iloc[:, -1].str.replace(';', '')
            df = df.loc[~df.index.duplicated(keep='first')].T
            dfs.append(df)
        return pd.concat(dfs).set_index(attr.index)

    def gtf_to_db(self, fpath, con):
        columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

        df_chunk = pd.read_csv(fpath, sep='\t', names=columns, comment='#')
        attr = df_chunk['attribute'].str.replace('"', '').str.split('; ')
        df_attr = self.get_attr(attr)
        df = pd.concat([df_chunk.drop('attribute', axis=1), df_attr], axis=1)

        for chr, df_chr in df.groupby('chromosome'):
            for str, df_str in df_chr.groupby('strand'):
                tname = '_'.join([chr, str])
                try:
                    df_str.to_sql(tname, con, if_exists='append', index=None)
                except:
                    df_org = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
                    df_str = pd.concat([df_str, df_org])
                    df_str.to_sql(tname, con, if_exists='replace', index=None)

    # def gtf_to_db(self, gtf_file, con):
    #     gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    #     dirname, fname = os.path.split(gtf_file)
    #
    #     df = pd.read_csv(gtf_file, names=gtf_columns, comment='#', sep='\t')
    #     df = df[df['feature'] == 'transcript'].reset_index(drop=True)
    #     df['chromosome'] = 'chr' + df['chromosome'].astype('str')
    #     attribute = df['attribute'].str.split('; ')
    #
    #     pkm = []
    #     for idx in attribute.index:
    #         attr = attribute[idx]
    #         dict = {}
    #         for a in attr:
    #             key, value = a.split(' ')
    #             dict[key] = value.replace('"', '')
    #
    #         row = [idx, None, None, None, dict['cov'], dict['FPKM'], dict['TPM'].replace(';', '')]
    #         if 'ref_gene_name' in dict:
    #             row[1] = dict['ref_gene_name']
    #             row[2] = dict['transcript_id']
    #             row[3] = dict['transcript_name']
    #         pkm.append(row)
    #
    #     df_res = pd.DataFrame(data=pkm, columns=['index', 'gene_name', 'transcript_id', 'transcript_name', 'cov', 'FPKM', 'TPM'])
    #     df_res = df_res.set_index('index', drop=True)
    #     df = pd.concat([df, df_res], axis=1)
    #
    #     df = df.drop(['attribute', 'frame'], axis=1)
    #     for chr, df_chr in df.groupby('chromosome'):
    #         for str, df_str in df_chr.groupby('strand'):
    #             tname = '_'.join([os.path.splitext(fname)[0], chr, str])
    #             df_str.sort_values(by=['start', 'end']).to_sql(tname, con, index=None)

    def correlation_replicates(self, df, cols):
        from itertools import combinations

        report = pd.Series()
        for comb in combinations(cols, 2):
            coeff = np.corrcoef(df[comb[0]].fillna(0), df[comb[1]].fillna(0))
            report['&'.join(comb)] = coeff[0, 1]
        return report

    def avg_rna_seq_by_tissues(self):
        df = pd.read_excel('RNA-seq_data_structure.xlsx')
        columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'gene_name', 'transcript_id', 'transcript_name']

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq.db')
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_tissue.db'))

        writer = pd.ExcelWriter(os.path.join(self.root, 'database/RNA-seq/out', 'repicates_duplicates.xlsx'), engine='xlsxwriter')
        for src, df_src in df.groupby('tissue'):
            dfs = []
            index = []
            for fid in df_src['fid']:
                df_fid = pd.read_sql_query("SELECT * FROM '{}' WHERE FPKM>''".format(fid), con)
                df_fid.index = df_fid['chromosome'] + ':' + df_fid['start'].astype(str) + '-' + df_fid['end'].astype(str)
                index.append(set(df_fid.index))
                dfs.append(df_fid)

            index = sorted(list(set.union(*index)))
            df_res = pd.DataFrame(index=index, columns=df_fid.columns)
            cols = []
            for i, df_fid in enumerate(dfs):
                df_fid['FPKM'] = df_fid['FPKM'].astype(float)
                cols.append('FPKM (Rep {})'.format(i + 1))
                df_res.loc[df_fid.index, cols[-1]] = df_fid['FPKM']
                df_res.loc[df_fid.index, columns] = df_fid[columns]

            # df_res = df_res.fillna(0)
            report = self.correlation_replicates(df_res, cols)
            report.to_excel(writer, sheet_name=src)

            df_res['FPKM'] = df_res[cols].mean(axis=1)
            df_res.to_sql(src, con_out, if_exists='replace', index=None)
        writer.save()
        writer.close()

    def avg_fantom_by_tissue(self):
        columns = ['chromosome', 'start', 'end', 'loc', 'score', 'strand']
        dirname_ = os.path.join(self.root, 'database/Fantom/v5/tissues')
        df = pd.read_excel('tissues_fantom_rna.xlsx', sheet_name='Sheet1')

        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'FANTOM_tissue.db')
        con_out = sqlite3.connect(fpath)

        writer = pd.ExcelWriter(os.path.join(self.root, 'database/Fantom/v5/tissues', 'repicates_duplicates.xlsx'), engine='xlsxwriter')
        for idx in df.index:
            tissue_rna = df.loc[idx, 'RNA-seq']
            tissue = df.loc[idx, 'FANTOM']
            dfs = []
            index = []

            dirname = os.path.join(dirname_, tissue)
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

            report = self.correlation_replicates(df_res, cols)
            report.to_excel(writer, sheet_name=tissue)
            df_res['score'] = df_res[cols].mean(axis=1)
            print(df_res[cols].sum(axis=0))
            df_res.drop(['loc'], axis=1).to_sql(tissue_rna, con_out, if_exists='replace', index=None)
        writer.save()
        writer.close()

    def avg_fantom_by_celllines(self):
        columns = ['chromosome', 'start', 'end', 'loc', 'score', 'strand']
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        df = pd.read_excel(os.path.join(dirname, 'list', '00_human.cell_line.hCAGE.hg19.assay_sdrf2.xlsx'))

        flist = os.listdir(dirname)
        flist = [x for x in flist if x.endswith('.gz')]
        flist_map = {}
        for fname in flist:
            fid = fname.split('.')[-5]
            flist_map[fid] = fname

        fpath = os.path.join(dirname, 'human_hCAGE.db')
        con_out = sqlite3.connect(fpath.replace('.db', '_celllines.db'))
        df_grp = df.groupby('cell line')

        df_report = pd.DataFrame(index=df_grp.groups, columns=['#data', 'fid'])
        writer = pd.ExcelWriter(os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'repicates_duplicates.xlsx'), engine='xlsxwriter')
        for src, df_src in df_grp:
            print(src)

            fid = ';'.join(df_src['Extract Name'])
            df_report.loc[src, '#data'] = df_src.shape[0]
            df_report.loc[src, 'fid'] = fid

            dfs = []
            index = []
            for idx in df_src.index:
                fid = df_src.loc[idx, 'Extract Name']
                fname = flist_map[fid]
                fpath = os.path.join(dirname, fname)

                df = pd.read_csv(fpath, sep='\t', names=columns, compression='gzip')
                df.index = df['chromosome'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
                index.append(set(df.index))
                dfs.append(df)

            if len(dfs) <= 1:
                dfs[0]['score'] = (dfs[0]['score'] / dfs[0]['score'].sum()) * 1e6
                dfs[0].to_sql(src, con_out, if_exists='replace', index=None)
            else:
                index = sorted(list(set.union(*index)))
                df_res = pd.DataFrame(index=index, columns=dfs[-1].columns)
                cols = []
                for i, df_fid in enumerate(dfs):
                    df_fid['score'] = df_fid['score'].astype(float)
                    df_fid['score'] = (df_fid['score'] / df_fid['score'].sum()) * 1e6
                    cols.append('score (Rep {})'.format(i + 1))
                    df_res.loc[df_fid.index, cols[-1]] = df_fid['score']
                    df_res.loc[df_fid.index, columns] = df_fid[columns]
                    print(df_res.loc[df_fid.index, columns].sum())

                df_res['score'] = df_res[cols].mean(axis=1)
                df_res.to_sql(src, con_out, if_exists='replace', index=None)

        writer.save()
        writer.close()
        df_report.to_excel('cell_line.xlsx')

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
    if con.hostname == 'mingyu-Precision-Tower-781':
        con.to_server()

    else:
        # con.stats_by_tissue()
        # con.avg_rna_seq_by_tissues()
        con.avg_fantom_by_tissue()

        # from Util import Util
        # from Database import Database
        #
        # ut = Util()
        # # fpath = os.path.join(ut.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        # fpath = os.path.join(ut.root, 'database/Fantom/v5/tissues', 'FANTOM_tissue.db')
        # tlist = Database.load_tableList(sqlite3.connect(fpath))
        # for tname in tlist:
        #     ut.split(fpath, tname)

        con.avg_fantom_by_celllines()

        # from Correlation_gpu import Correlation
        # cor = Correlation()
        # cor.correlation_fan_rna()
        # cor.high_correlation()

        # fpath = os.path.join(con.root, 'database/gencode', 'gencode.v32lift37.annotation.gtf')
        # conn = sqlite3.connect(fpath.replace('.gtf', '.db'))
        # con.gtf_to_db(fpath, conn)

        # dirname = os.path.join(con.root, 'database', 'RNA-seq')
        # con.gtf_to_db_all(dirname)

        # con.avg_rna_seq_by_celllines()
        # dirname = os.path.join(con.root, 'database', 'RNA-seq', 'bam')
        # flist = os.listdir(dirname)
        # for fname in flist:
        #     fpath = os.path.join(dirname, fname)
        #     con.bam_to_gtf(fpath)
        #     con.gtf_to_db(fpath.replace('.bam', '.gtf'))
