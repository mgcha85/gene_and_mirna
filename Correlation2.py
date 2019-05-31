import pandas as pd
import numpy as np
import sqlite3
import os
import socket
import sys
from scipy.stats import spearmanr

from Database import Database
from Server import Server


class Correlation2:
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
        self.cells = []

    def to_server(self):
        server = Server()
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server.job_script(fname, time='03:00:00')

        server_root = os.path.join(server.server, 'source/gene_and_mirna')
        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def get_reference(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql_query("SELECT * FROM 'human_genes'", con)
        df_unique = df.drop_duplicates(subset=['gene_name'], keep=False)

        idx_duplicates = list(set(df.index) - set(df_unique.index))
        df_duplicates = df.loc[idx_duplicates]
        df_duplicates = df_duplicates[df_duplicates['source'] == 'ENSEMBL']
        df_duplicates = df_duplicates.drop_duplicates(subset=['gene_name'], keep=False)
        return pd.concat([df_unique, df_duplicates]).reset_index(drop=True)

    def pre_data(self):
        fpath = os.path.join(self.root, 'software/stringtie-1.3.6', 'genes.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        out_con = sqlite3.connect(fpath.replace('.db', '_2.db'))

        for tname in tlist:
            df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
            attribute = df_rna['attribute'].str.split('; ')

            pkm = []
            for attr in attribute:
                dict = {}
                for a in attr[:-1]:
                    key, value = a.split(' ')
                    dict[key] = value.replace('"', '')

                if 'gene_name' not in dict:
                    pkm.append(np.nan)
                    continue
                pkm.append(dict['gene_name'])

            df_res = pd.DataFrame(data=pkm, columns=['gene_name'])
            df_rna = pd.concat([df_rna, df_res], axis=1)
            df_rna.drop(['attribute', 'frame'], axis=1).to_sql(tname, out_con, index=None, if_exists='append')

    def fantom_unique_gene(self):
        fpath_rna = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue_out.db')
        con_rna = sqlite3.connect(fpath_rna)
        con_out = sqlite3.connect(fpath_rna.replace('.db', '2.db'))
        tlist_rna = Database.load_tableList(con_rna)

        for tname in tlist_rna:
            print(tname)
            df_rna = pd.read_sql_query("SELECT * FROM {}".format(tname), con_rna)
            tot = df_rna['score'].sum()
            df_grp = df_rna.groupby('gene_name')

            contents = []
            for gene, df_sub in df_grp:
                score = df_sub['score'].sum() / tot
                chromosome = df_sub['chromosome'].iloc[0]
                strand = df_sub['strand'].iloc[0]
                tss = ';'.join(df_sub['start'].astype(str))
                scores = ';'.join(df_sub['score'].astype(str))
                contents.append([chromosome, strand, tss + ',' + scores, score, gene])

            df_res = pd.DataFrame(contents, columns=['chromosome', 'strand', 'tss', 'score', 'gene_name'])
            df_res.to_sql(tname, con_out, if_exists='replace', index=None)

    def rna_unique_gene(self):
        fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_out.db')
        con_rna = sqlite3.connect(fpath_rna)
        con_out = sqlite3.connect(fpath_rna.replace('.db', '2.db'))
        tlist_rna = Database.load_tableList(con_rna)

        for tname in tlist_rna:
            print(tname)
            df_rna = pd.read_sql_query("SELECT * FROM {} WHERE FPKM>0".format(tname), con_rna)
            df_rna['FPKM'] = df_rna['FPKM'].astype(float)
            df_grp = df_rna.groupby('gene_name')

            contents = []
            for gene, df_sub in df_grp:
                midx = df_sub['FPKM'].idxmax()
                contents.append(df_sub.loc[midx:midx, :])
            df_res = pd.concat(contents)
            df_res.to_sql(tname, con_out, if_exists='replace', index=None)

    def get_valid_tissues(self, tlist_rna, tlist_fan):
        table_names = {}
        for tname_fan in tlist_fan:
            tlist_rna_sub = [x for x in tlist_rna if tname_fan in x]
            table_names[tname_fan] = tlist_rna_sub
        return table_names

    def correlation(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'correlation.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        dfs = {}
        gene_names = []
        for tname in tlist:
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con, index_col='gene_name')
            if df.empty:
                continue
            dfs[tname] = df
            gene_names.append(set(df.index))

        report = []
        gene_names = sorted(list(set.intersection(*gene_names)))
        for i, gname in enumerate(gene_names):
            if i % 100 == 0:
                print('{:0.2f}%'.format(100 * (i + 1) / len(gene_names)))

            scores = []
            fpkms = []
            tissues = []
            tsss = []
            loci = ''
            tissue = ''
            for j, (tname, df) in enumerate(dfs.items()):
                tissue__ = tname.split('_')[0]
                if tissue != tissue__:
                    tissue = tissue__
                    tissues.append(tissue)
                    scores.append(df.loc[gname, 'score'])
                    fpkms.append(df.loc[gname, 'FPKM'])
                    tsss.append(df.loc[gname, 'tss'])
                if j == 0:
                    chromosome = df.loc[gname, 'chromosome']
                    strand = df.loc[gname, 'strand']
                    if strand == '+':
                        tss = df.loc[gname, 'start']
                    else:
                        tss = df.loc[gname, 'end']
                    loci = ':'.join([chromosome, str(tss), strand])

            if i == 0:
                print(tissues)
            fpkm_str = [str(x) for x in fpkms]
            spearman_xcor = spearmanr(scores, fpkms)
            report.append([gname, loci, ':'.join(tsss), ','.join(fpkm_str), spearman_xcor.correlation])

        df_rep = pd.DataFrame(data=report, columns=['gnames', 'tss', 'scores', 'fpkms', 'corr (xcor)'])
        df_rep.to_excel('correlation.xlsx', index=None)

    def run(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con = sqlite3.connect(fpath)
        df_ref = pd.read_sql_query("SELECT chromosome, start, end, strand, gene_name FROM 'human_genes_wo_duplicates'", con, index_col='gene_name')

        fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue_out2.db')
        con_rna = sqlite3.connect(fpath_rna)
        tlist_rna = Database.load_tableList(con_rna)

        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/tissues/out/', 'fantom_cage_by_tissue_out2.db')
        con_fan = sqlite3.connect(fpath_fan)
        tlist_fan = Database.load_tableList(con_fan)

        out_path = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'correlation.db')
        out_con = sqlite3.connect(out_path)

        tlist = self.get_valid_tissues(tlist_rna, tlist_fan)
        for tissue, tname_rna in tlist.items():
            print(tissue)
            df_fan = pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con_fan, index_col='gene_name')
            for rtname in tname_rna:
                df_rna = pd.read_sql_query("SELECT * FROM '{}'".format(rtname), con_rna, index_col='gene_name')

                common_genes = sorted(list(set.intersection(set(df_fan.index), set(df_rna.index))))
                df_res = pd.concat([df_ref.loc[common_genes], df_fan.loc[common_genes, 'score'], df_rna.loc[common_genes, 'FPKM'], df_fan.loc[common_genes, 'tss']], axis=1)
                df_res.to_sql(rtname, out_con, if_exists='replace')


if __name__ == '__main__':
    from Comparison_gpu import Comparison
    comp = Comparison()

    cor = Correlation2()
    if cor.hostname == 'mingyu-Precision-Tower-7810':
        # cor.correlation()
        cor.to_server()

    else:
        comp.fantom_to_gene(500)
        cor.fantom_unique_gene()
        # cor.rna_unique_gene()
        cor.run()
        cor.correlation()
