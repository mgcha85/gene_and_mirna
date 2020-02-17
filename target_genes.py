import sqlite3
import os
import pandas as pd
import socket
from Database import Database


class target_genes:
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

        which = 'stokes'
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

    def group_by_mir(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed.db')
        con = sqlite3.connect(fpath)
        # tlist = Database.load_tableList(con)

        for tname in ['miranda', 'rna22', 'ts']:
            print(tname)
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            contents = []
            for mir, df_sub in df.groupby('miRNA'):
                genes = ';'.join(sorted(list(set(df_sub['genes']))))
                contents.append([mir, genes])
            df_res = pd.DataFrame(contents, columns=['mir', 'genes'])
            df_res.to_sql(tname + '_grp_mir', con, index=None, if_exists='replace')

    def union(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed.db')
        con = sqlite3.connect(fpath)
        tlist = [x for x in Database.load_tableList(con) if '_grp_mir' in x]

        mir = set()
        dfs = []
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='mir')
            df = df.loc[~df.index.duplicated(keep='first')]
            mir = set.union(mir, set(df.index))
            dfs.append(df)

        df_res = pd.DataFrame(index=sorted(list(mir)), columns=range(len(dfs)))
        for i, df in enumerate(dfs):
            df_res.loc[df.index, i] = df.loc[df.index, 'genes']

        result = []
        for mir in df_res.index:
            union = []
            for col in df_res.columns:
                if isinstance(df_res.loc[mir, col], str):
                    genes = set(df_res.loc[mir, col].split(';'))
                    union.append(genes)
            union = set.union(*union)
            result.append([mir, ';'.join(sorted(list(union)))])
        df_res = pd.DataFrame(result, columns=['mir', 'genes'])
        df_res.to_sql('union', con, index=None, if_exists='replace')

    def intersection(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed.db')
        con = sqlite3.connect(fpath)
        tlist = [x for x in Database.load_tableList(con) if '_grp_mir' in x]

        mir = []
        dfs = []
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='mir')
            df = df.loc[~df.index.duplicated(keep='first')]
            mir.append(set(df.index))
            dfs.append(df)

        mir = set.intersection(*mir)
        df_res = pd.DataFrame(index=sorted(list(mir)), columns=range(len(dfs)))
        for i, df in enumerate(dfs):
            df_res.loc[mir, i] = df.loc[mir, 'genes']

        result = []
        for mir in df_res.index:
            inter = []
            for col in df_res.columns:
                inter.append(set(df_res.loc[mir, col].split(';')))
            inter = set.intersection(*inter)
            result.append([mir, ';'.join(sorted(list(inter)))])
        df_res = pd.DataFrame(result, columns=['mir', 'genes'])
        df_res.to_sql('intersection', con, index=None, if_exists='replace')

    def make_table(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed.db')
        con = sqlite3.connect(fpath)

        fpath_ref = os.path.join(self.root, 'database/gencode', 'gencode.v32lift37.annotation_attr_merged.db')
        con_ref = sqlite3.connect(fpath_ref)
        df_ref = pd.read_sql("SELECT chromosome, start, end, strand, transcript_id, gene_name FROM 'gencode.v32lift37' WHERE feature='transcript' AND "
                "gene_type='protein_coding' AND transcript_type='protein_coding'", con_ref)
        df_ref_grp = df_ref.groupby('gene_name')
        df = pd.read_sql("SELECT * FROM 'intersection'", con, index_col='mir')

        contents = {}
        for i, mir in enumerate(df.index):
            if (i + 1) % 100 == 0 or i + 1 == df.shape[0]:
                print('{} / {}'.format(i+1, df.shape[0]))
            for gene in df.loc[mir, 'genes'].split(';'):
                if gene in contents or gene not in df_ref_grp.groups:
                    continue
                df_gene = df_ref_grp.get_group(gene)
                contents[gene] = df_gene.loc[df_gene.index[0]: df_gene.index[0]]

        fpath_out = os.path.join(self.root, 'database/gencode', 'other_researches_inter_fan_rna_100.db')
        con_out = sqlite3.connect(fpath_out)
        for gene, row in contents.items():
            row.drop(['chromosome', 'strand'], axis=1).to_sql('_'.join([row['chromosome'].iloc[0], row['strand'].iloc[0]]), con_out,
                                                              index=None, if_exists='append')


if __name__ == '__main__':
    tg = target_genes()
    if tg.hostname == 'mingyu-Precision-Tower-7810':
        tg.to_server()
    else:
        # tg.group_by_mir()
        # tg.union()
        # tg.intersection()
        tg.make_table()
