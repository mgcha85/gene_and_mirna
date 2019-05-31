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
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        self.bed_columns = ['chromosome', 'start', 'end', 'location', 'score', 'strand']

    def to_server(self):
        server = Server()
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server.job_script(fname, time='00:30:00')

        server_root = os.path.join(server.server, 'source/gene_and_mirna')
        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def gtf_to_db(self):
        # dirname = os.path.join(self.root, 'database/RNA-seq/gtf')
        # flist = [os.path.join(dirname, x) for x in os.listdir(dirname)]

        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        flist = []
        for path, subdirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.gtf'):
                    flist.append(os.path.join(path, name))

        for fpath in flist:
            dirname, fname = os.path.split(fpath)
            df = pd.read_csv(fpath, sep='\t', names=self.gtf_columns, comment='#')
            df = df[df['feature'] == 'transcript']
            df = df.reset_index(drop=True)
            # df['chromosome'] = 'chr' + df['chromosome'].astype('str')
            attribute = df['attribute'].str.split('; ')

            pkm = []
            for attr in attribute:
                dict = {}
                for a in attr:
                    key, value = a.split(' ')
                    dict[key] = value.replace('"', '')
                pkm.append([dict['FPKM'], dict['TPM'].replace(';', '')])

            fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue_out.db')
            con = sqlite3.connect(fpath)
            df_res = pd.DataFrame(data=pkm, columns=['FPKM', 'TPM'])
            df = pd.concat([df, df_res], axis=1)
            tname = os.path.splitext(fname)[0].split('%')[0]
            df.drop(['frame', 'attribute'], axis=1).to_sql(tname, con, index=None, if_exists='replace')

    def bed_to_db(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        flist = []
        for path, subdirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.gz'):
                    flist.append(os.path.join(path, name))

        out_path = os.path.join(dirname, 'out', 'fantom_cage_by_tissue.db')
        con = sqlite3.connect(out_path)
        for fpath in flist:
            df = pd.read_csv(fpath, sep='\t', names=self.bed_columns, comment='#', compression='gzip')

            dirname, fname = os.path.split(fpath)
            tissue = dirname.split('/')[-1]
            df.to_sql(tissue, con, index=None, if_exists='append')

    def get_merge_table(self, fpath):
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        dfs = []
        for tname in tlist:
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            dfs.append(df)
        return pd.concat(dfs)

    def get_genes(self, fpath):
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        addition = np.zeros(2)

        for tname in tlist:
            # df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            df = pd.read_sql_query("SELECT * FROM '{}' WHERE fan_start>0 AND fan_end>0".format(tname), con)
            gene_names = np.unique(df['gene_name'].values)
            transcript_name = np.unique(df['transcript_name'].values)
            addition[0] += gene_names.shape[0]
            addition[1] += transcript_name.shape[0]
        return addition

    def get_fantom(self, fpath):
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        addition = 0

        for tname in tlist:
            df = pd.read_sql_query("SELECT * FROM '{}' WHERE fan_start>0 AND fan_end>0".format(tname), con)
            addition += df.shape[0]
        return addition

    def protein_coding(self):
        fpaths = {'GEN': os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation.db'),
                  'ENS': os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19.db'),
                  'UCSC': os.path.join(self.root, 'database/UCSC/Genes', 'genes.db')}

        out = np.zeros((len(fpaths), 2))
        df = pd.DataFrame(data=out, index=fpaths.keys(), columns=['org', 'protein_coding'])
        for i, (key, fpath) in enumerate(fpaths.items()):
            dirname, fname = os.path.split(fpath)

            df.iloc[i]['org'] = self.get_merge_table(fpath).shape[0]
            fname_, ext = os.path.splitext(fname)
            df.iloc[i]['protein_coding'] = self.get_merge_table(os.path.join(dirname, fname_ + '_pc' + ext)).shape[0]

        print(df)

    def gene_name(self):
        fpaths = {'GEN': os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation_pc.db'),
                  'ENS': os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19_pc.db'),
                  'UCSC': os.path.join(self.root, 'database/UCSC/Genes', 'genes_pc.db')}

        out = np.zeros((len(fpaths), 2))
        df = pd.DataFrame(data=out, index=fpaths.keys(), columns=['#gene_names', '#transctipt_names'])
        for i, (key, fpath) in enumerate(fpaths.items()):
            df.iloc[i, :] = self.get_genes(fpath)
        print(df)

    def after_fantom(self):
        fpaths = {'GEN': os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation_pc.db'),
                  'ENS': os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19_pc.db'),
                  'UCSC': os.path.join(self.root, 'database/UCSC/Genes', 'genes_pc.db')}

        out = np.zeros((len(fpaths), 1))
        df = pd.DataFrame(data=out, index=fpaths.keys(), columns=['#after_fantom'])
        for i, (key, fpath) in enumerate(fpaths.items()):
            df.iloc[i, :] = self.get_fantom(fpath)
        print(df)

    def get_cage_peak_by_id(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        df_tss = pd.read_sql_query("SELECT * FROM 'TSS_human'", con)
        df_gene = pd.read_sql_query("SELECT * FROM 'human_gene'", con, index_col='gene')
        
        df_tss.loc[:, 'promoter'] = df_tss['promoter'].str.split(',', expand=True)[0]
        df_tss = df_tss.set_index('promoter', drop=True)

        root_dir = os.path.join(self.root, 'database/Fantom/v5/tissues')
        dir_list = os.listdir(root_dir)
        dir_list = [x for x in dir_list if os.path.isdir(os.path.join(root_dir, x))]

        # load cage peak id by tissue
        df_tissue = {}
        for tissue in dir_list:
            dirname = os.path.join(root_dir, tissue)
            flist = [x for x in os.listdir(dirname) if x.endswith('.csv')]
            for fname in flist:
                if '~' in fname:
                    continue

                fpath = os.path.join(dirname, fname)
                df = pd.read_csv(fpath, index_col=0)
                N = len(flist)
                if tissue in df_tissue:
                    df = df[df['TPM'] > 1]
                    int_idx = set.intersection(set(df_tissue[tissue].index), set(df.index))
                    df_tissue[tissue].loc[int_idx, 'TPM'] += df.loc[int_idx, 'TPM'] / N
                else:
                    df_tissue[tissue] = df[df['TPM'] > 1] / N

        # get gene locations
        not_in_list = []
        out_con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5', 'hg19.tissue.db'))

        for tissue, df in df_tissue.items():
            print(tissue, df.shape[0])
            result = []

            for cage_peak in df.index:
                if cage_peak not in df_tss.index:
                    not_in_list.append(cage_peak)
                else:
                    gene_name = df_tss.loc[cage_peak, 'gene']
                    tpm = df.loc[cage_peak, 'TPM']
                    if gene_name in df_gene.index:
                        gene_loc = df_gene.loc[gene_name, ['start', 'end']]
                        result.append([*df_tss.loc[cage_peak, ['chromosome', 'start', 'end', 'strand', 'gene']],
                                       *gene_loc, cage_peak, tpm])
                    else:
                        print(gene_name)

            df_res = pd.DataFrame(data=result, columns=['chromosome', 'start', 'end', 'strand', 'gene_name',
                                                        'gene_start', 'gene_end', 'cage_peak', 'TPM'])
            df_res.sort_values(by=['chromosome', 'start']).to_sql(tissue, out_con, if_exists='replace', index=None)

        with open('not_in_list.txt', 'wt') as f:
            f.write('\n'.join(not_in_list))


if __name__ == '__main__':
    dp = data_preparation()
    if dp.hostname == 'mingyu-Precision-Tower-7810':
        dp.to_server()
    else:
        dp.bed_to_db()
    # dp.gtf_to_db()
