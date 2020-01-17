import pandas as pd
import numpy as np
import os
import sqlite3
import socket
from Database import Database


class miRNA:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
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

        server.job_script(fname, src_root=server_root, time='08:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def convert(self):
        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.bed")
        df = pd.read_csv(fpath, sep='\t', index_col=0, names=['locations'])
        df.index.name = 'miRNA'

        out_path = fpath.replace('.bed', '.db')
        con = sqlite3.connect(out_path)
        df.to_sql('original', con)

    def seperate_loc(self):
        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.db")
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'original'", con)
        df['locations'] = df['locations'].str.split(';')

        df_res = pd.DataFrame(index=df.index, columns=['miRNA', 'chromosome', 'starts', 'ends'])
        for idx in df.index:
            mir = df.loc[idx, 'miRNA']
            loc = df.loc[idx, 'locations']
            N = len(loc)
            df_each = pd.DataFrame(data=np.zeros((N, 3)), columns=['chromosome', 'start', 'end'], index=range(N))

            for i, l in enumerate(loc):
                chr, start, end = l.split('_')
                df_each.loc[i, 'chromosome'] = chr
                df_each.loc[i, 'start'] = int(start)
                df_each.loc[i, 'end'] = int(end)

            df_res.loc[idx, 'miRNA'] = mir
            df_res.loc[idx, 'chromosome'] = chr
            df_res.loc[idx, 'starts'] = ';'.join(df_each['start'].astype(int).astype(str))
            df_res.loc[idx, 'ends'] = ';'.join(df_each['end'].astype(int).astype(str))
        df_res.to_sql('seperate_tss', con, index=None)

    def mean_tss(self):
        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.db")
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'seperate_tss'", con, index_col='miRNA')
        df = df.loc[~df.index.duplicated(keep='first')]

        df_res = pd.DataFrame(index=df.index, columns=['chromosome', 'tss'])
        for idx in df.index:
            starts = df.loc[idx, 'starts'].split(';')
            ends = df.loc[idx, 'ends'].split(';')

            smin = min(map(int, starts))
            emin = min(map(int, ends))
            smax = max(map(int, starts))
            emax = max(map(int, ends))

            tss = (min([smin, emin]) + max([smax, emax])) // 2
            df_res.loc[idx, 'chromosome'] = df.loc[idx, 'chromosome']
            df_res.loc[idx, 'tss'] = tss
        df_res.to_sql('mean_tss', con)

    def add_strand(self):
        con_ref = sqlite3.connect(os.path.join(self.root, 'database', 'papers.db'))
        tlist = Database.load_tableList(con_ref)
        tlist = [x for x in tlist if x != 'miRNA_Type']

        fpath = os.path.join(self.root, "database", "consistent_miRNA_330.db")
        con = sqlite3.connect(fpath)

        df = pd.read_sql("SELECT * FROM 'mean_tss'", con)
        for idx in df.index:
            mirname = df.loc[idx, 'miRNA']
            print(mirname)

            for tname in tlist:
                columns = Database.load_table_columns(con_ref, tname)
                try:
                    sidx = columns.index('strand')
                except:
                    sidx = columns.index('Strand')

                df_ref = pd.read_sql("SELECT miRNA, {} FROM '{}' WHERE miRNA='{}'".format(columns[sidx], tname, mirname), con_ref, index_col='miRNA')
                if not df_ref.empty:
                    if df_ref.shape[0] > 1:
                        df_ref = df_ref.iloc[:1]
                    strand = df_ref.loc[mirname, columns[sidx]]
                    df.loc[idx, 'strand'] = strand
                    break
        df.to_sql('mean_tss', con, if_exists='replace', index=None)

    def from_rna_seq(self):
        from Database import Database

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con = sqlite3.connect(fpath)

        dfs = []
        for tname in Database.load_tableList(con):
            df = pd.read_sql("SELECT * FROM '{}' WHERE ref_gene_name LIKE '%MIR%'".format(tname), con, index_col='reference_id')
            df['tissue'] = tname
            dfs.append(df)
        df = pd.concat(dfs)
        df = df.drop_duplicates('ref_gene_name')

        out_path = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_mir.db')
        out_con = sqlite3.connect(out_path)
        df.to_sql('MIR', out_con)

    def convert_name(self):
        fpath = os.path.join(self.root, 'database/mirbase/22', 'mirbase.db')
        con = sqlite3.connect(fpath)
        df_mir = pd.read_sql("SELECT * FROM 'hsa_liftover_hg19' WHERE feature='miRNA_primary_transcript'", con)

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_mir.db')
        con = sqlite3.connect(fpath)
        df_rna = pd.read_sql("SELECT * FROM 'MIR'", con, index_col='reference_id')

        # 1 by name
        idx = df_rna[df_rna['ref_gene_name'].str.contains('hsa-mir')].index
        df_rna.loc[idx, 'Name'] = df_rna.loc[idx, 'ref_gene_name']

        # 2 by position
        df_mir.index = df_mir['chromosome'] + ':' + df_mir['start'].astype(str) + '-' + df_mir['end'].astype(str) + ',' + df_mir['strand']
        df_rna.index = df_rna['chromosome'] + ':' + df_rna['start'].astype(str) + '-' + df_rna['end'].astype(str) + ',' + df_rna['strand']

        df_mir.index = df_mir['chromosome'] + ':' + df_mir['start'].astype(str) + '-' + df_mir['end'].astype(str) + ',' + df_mir['strand']
        df_rna.index = df_rna['chromosome'] + ':' + df_rna['start'].astype(str) + '-' + df_rna['end'].astype(str) + ',' + df_rna['strand']

        inter = set.intersection(set(df_mir.index), set(df_rna.index))

        df_mir = df_mir.loc[~df_mir.index.duplicated(keep='first')]
        df_rna = df_rna.loc[~df_rna.index.duplicated(keep='first')]
        df_rna.loc[inter, 'Name'] = df_mir.loc[inter, 'Name']
        df_rna_null = df_rna[df_rna['Name'].isnull()]

        # 3 by different name
        df_mir['mirName'] = df_mir['Name'].str.replace('hsa-', '').str.replace('-', '').str.replace('mir', '')
        df_mir = df_mir.set_index('mirName', drop=True)

        df_rna_null['mirName'] = df_rna_null['ref_gene_name'].str.replace('-', '').str.replace('MIR', '').str.lower()
        df_rna_null['mirName'] = df_rna_null['mirName'].str.replace('hg', '')
        df_rna_null = df_rna_null.set_index('mirName', drop=True)

        inter = set.intersection(set(df_mir.index), set(df_rna_null.index))
        inter = sorted(list(inter))
        df_mir = df_mir.loc[~df_mir.index.duplicated(keep='first')]

        df_rna_null.loc[inter, 'Name'] = df_mir.loc[inter, 'Name']
        df_rna_null.index = df_rna_null['chromosome'] + ':' + df_rna_null['start'].astype(str) + '-' + df_rna_null['end'].astype(str) + ',' + df_rna_null['strand']
        df_rna.loc[df_rna_null.index, 'Name'] = df_rna_null.loc[:, 'Name']
        df_rna.sort_values(by=['chromosome', 'start']).to_sql('MIR_NAME', con, index=None, if_exists='replace')

    def get_consistent_mir(self):
        fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_mir.db')
        con_rna = sqlite3.connect(fpath_rna)
        df_rna = pd.read_sql("SELECT * FROM 'MIR_NAME'", con_rna, index_col='Name')

        fpath = os.path.join(self.root, 'database', 'consistent_miRNA_330.db')
        con = sqlite3.connect(fpath)
        df_con = pd.read_sql("SELECT * FROM 'mean_tss'", con, index_col='miRNA')

        common = set.intersection(set(df_rna.index), set(df_con.index))
        df_rna.loc[common].to_sql('MIR_consistent', con_rna, if_exists='replace')

    def get_expression(self):
        fpath_rna = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_mir.db')
        con_rna = sqlite3.connect(fpath_rna)
        df_rna = pd.read_sql("SELECT * FROM 'MIR_consistent' GROUP BY Name HAVING MAX(FPKM)", con_rna, index_col='Name')

        fpath_tis = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con_tis = sqlite3.connect(fpath_tis)

        tissues = Database.load_tableList(con_tis)
        df_buffer = pd.DataFrame(np.zeros((df_rna.shape[0], len(tissues))), index=df_rna.index, columns=tissues)
        for i, mir in enumerate(df_rna.index):
            if i % 20 == 0 or i + 1 == df_rna.shape[0]:
                print('{} / {}'.format(i+1, df_rna.shape[0]))

            rid = df_rna.loc[mir, 'ref_gene_id']
            for tis in tissues:
                sql = "SELECT ref_gene_id, FPKM FROM (SELECT ref_gene_id, FPKM FROM '{}' WHERE ref_gene_id='{}') GROUP BY ref_gene_id HAVING MAX(FPKM)".format(tis, rid)
                df_tis = pd.read_sql(sql, con_tis, index_col='ref_gene_id')
                if df_tis.empty:
                    continue
                df_buffer.loc[mir, tis] = df_tis.loc[rid, 'FPKM']
        df_buffer.to_sql('MIR_expression', con_rna, if_exists='replace')

    def get_gene_expression(self):
        fpath_gen = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_100.db')
        con_gen = sqlite3.connect(fpath_gen)
        dfs = []
        for tname in Database.load_tableList(con_gen):
            dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con_gen, index_col='transcript_id'))
        df_gen = pd.concat(dfs)

        fpath_tis = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con_tis = sqlite3.connect(fpath_tis)

        fpath_out = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_gene.db')
        con_out = sqlite3.connect(fpath_out)

        tissues = Database.load_tableList(con_tis)
        df_buffer = pd.DataFrame(np.zeros((df_gen.shape[0], len(tissues))), index=df_gen.index, columns=tissues)
        for i, tid in enumerate(df_gen.index):
            if i % 20 == 0 or i + 1 == df_gen.shape[0]:
                print('{} / {}'.format(i+1, df_gen.shape[0]))

            for tis in tissues:
                sql = "SELECT reference_id, FPKM FROM (SELECT reference_id, FPKM FROM '{}' WHERE reference_id='{}') GROUP BY reference_id HAVING MAX(FPKM)".format(tis, tid)
                df_tis = pd.read_sql(sql, con_tis, index_col='reference_id')
                if df_tis.empty:
                    continue
                df_buffer.loc[tid, tis] = df_tis.loc[tid, 'FPKM']
        df_buffer.to_sql('gene_expression', con_out, if_exists='replace')


if __name__ == '__main__':
    mir = miRNA()
    if mir.hostname == 'mingyu-Precision-Tower-781':
        mir.to_server()
    else:
        # mir.convert()
        # mir.seperate_loc()
        # mir.mean_tss()
        mir.get_gene_expression()
