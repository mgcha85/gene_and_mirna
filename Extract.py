import sqlite3
import os
import pandas as pd
import numpy as np
import socket
import sys
import multiprocessing
from Server import Server


class Extract:
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
        self.num_cores = multiprocessing.cpu_count()
        self.scope = 500

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

    def processInput(self, df, i, mir):
        if i % 1000 == 0 or i + 1 == df.shape[0]:
            print('{:0.2f}%'.format(100 * (i + 1) / df.shape[0]))

        mir_path = os.path.join(self.root, 'database', 'fantom5.db')
        con_mir = sqlite3.connect(mir_path)

        gene_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con_gene = sqlite3.connect(gene_path)

        mir_low = mir.lower()
        df_mir = pd.read_sql("SELECT * FROM 'human_promoters_wo_duplicates' WHERE premiRNA='{}'".format(mir_low),
                                   con_mir, index_col='premiRNA')
        if df_mir.empty:
            return

        mir_chr = df_mir.loc[mir_low, 'chromosome']
        mir_strand = df_mir.loc[mir_low, 'strand']

        gene = df.loc[mir, 'Target Gene']
        distance = None
        df_gene = pd.read_sql(
            "SELECT * FROM 'human_genes_wo_duplicates' WHERE gene_name='{}'".format(gene), con_gene,
            index_col='gene_name')
        if df_gene.empty:
            return

        gene_chr = df_gene.loc[gene, 'chromosome']
        gene_strand = df_gene.loc[gene, 'strand']
        if gene_strand == '+':
            offset = 'start'
        else:
            offset = 'end'

        if gene_strand == mir_strand and gene_chr == mir_chr:
            distance = abs(df_mir.loc[mir, offset] - df_gene.loc[gene, offset])

        return [mir, gene, mir_chr, mir_strand, gene_chr, gene_strand, distance]

    def get_distance(self):
        # from joblib import Parallel, delayed
        # import multiprocessing
        # num_cores = multiprocessing.cpu_count()
        # if num_cores > 8:
        #     num_cores = 8

        mir_path = os.path.join(self.root, 'database', 'fantom5.db')
        con_mir = sqlite3.connect(mir_path)

        gene_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con_gene = sqlite3.connect(gene_path)

        fpath = os.path.join(self.root, 'database', 'miRTarBase.db')
        con = sqlite3.connect(fpath)

        df = pd.read_sql("SELECT * FROM 'miRTarBase' WHERE Species_miRNA='Homo sapiens'", con)
        df = df.drop_duplicates(subset=['miRNA', 'Target Gene'])
        df = df.set_index('miRNA', drop=True)

        contents = []
        for i, mir in enumerate(df.index):
            if i % 1000 == 0 or i + 1 == df.shape[0]:
                print('{:0.2f}%'.format(100 * (i + 1) / df.shape[0]))

            mir_low = mir.lower()
            df_mir = pd.read_sql("SELECT * FROM 'human_promoters_wo_duplicates' WHERE miRNA='{}'".format(mir_low), con_mir, index_col='miRNA')
            df_mir = df_mir.loc[~df_mir.index.duplicated(keep='first')]
            if df_mir.empty:
                continue

            mir_chr = df_mir.loc[mir_low, 'chromosome']
            mir_strand = df_mir.loc[mir_low, 'strand']

            genes = df.loc[mir, 'Target Gene']
            for gene in genes:
                distance = None
                df_gene = pd.read_sql(
                    "SELECT * FROM 'human_genes_wo_duplicates' WHERE gene_name='{}'".format(gene), con_gene,
                    index_col='gene_name')
                if df_gene.empty:
                    continue

                gene_chr = df_gene.loc[gene, 'chromosome']
                gene_strand = df_gene.loc[gene, 'strand']
                if gene_strand == '+':
                    offset = 'start'
                else:
                    offset = 'end'

                if gene_strand == mir_strand and gene_chr == mir_chr:
                    distance = abs(df_mir.loc[mir, offset] - df_gene.loc[gene, offset])

                contents.append([mir, gene, mir_chr, mir_strand, gene_chr, gene_strand, distance])

        df_res = pd.DataFrame(contents, columns=['mir', 'gene', 'mir_chr', 'mir_strand', 'gene_chr', 'gene_strand', 'distance'])
        # df_res = df_res.dropna(how='all', axis=0)

        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector_corr3.xlsx'.format(self.scope))
        df_res.to_excel(out_path, index=None)

    def stats(self):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector_corr3.xlsx'.format(self.scope))
        df = pd.read_excel(fpath)
        df_sub = df[df['distance']]
        print('{:0.2f}%'.format(100 * df_sub.shape[0] / df.shape[0]))
        df_sub['distance'].hist(bins=100)
        plt.show()

    def get_high_corr(self, df):
        contents = []
        for mir in df.index:
            row = df.loc[mir, :]
            row = row[row > 0.6]
            if len(row) > 0:
                contents.append([mir, ';'.join(row.index), ';'.join(row.round(2).astype(str))])
        return pd.DataFrame(data=contents, columns=['miRNA', 'GENEs', 'corr (pearson)'])

    def split_mir_table(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        out_con = sqlite3.connect(fpath.replace('.db', '_mir.db'))

        df = pd.read_sql("SELECT * FROM 'human_promoters_wo_duplicates'", con)
        print(df.shape[0])

        df_chr = df.groupby('chromosome')
        cnt = 0
        for chr, df_sub in df_chr:
            df_sub_str = df_sub.groupby('strand')
            for str, df_sub_sub in df_sub_str:
                cnt += df_sub_sub.shape[0]
                print(cnt)
                df_sub_sub = df_sub_sub.rename(columns={"premiRNA": "gene_name"})
                if str == '+':
                    df_sub_sub['start'] = df_sub_sub['tss']
                else:
                    df_sub_sub['end'] = df_sub_sub['tss']
                df_sub_sub.to_sql('human_promoters_{}_{}'.format(chr, str), out_con, index=None, if_exists='replace')

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_corr.xlsx'.format(self.scope))
        df = pd.read_excel(fpath, index_col=0).T
        df = df.dropna(axis=1)
        df = self.get_high_corr(df)
        df.to_excel(fpath.replace('.xlsx', '2.xlsx'), index=None)


if __name__ == '__main__':
    ex = Extract()
    if ex.hostname == 'mingyu-Precision-Tower-7810':
        # ex.get_distance()
        ex.run()
    else:
        ex.run()
