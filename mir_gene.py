import sqlite3
import pandas as pd
import os
import socket
from scipy.stats import hypergeom
import numpy as np
import math
from Database import Database


class mir_gene:
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

    def to_sql(self):
        fpath = "/home/mingyu/Bioinformatics/database/feature_info.txt"
        df = pd.read_csv(fpath, sep='\t', names=['Target Gene', 'miRNA', 'binding energy', 'structural accessibility', 'context score'])
        df['Target Gene'] = df['Target Gene'].str.split('|', expand=True)[0]
        df['Target Gene'] = df['Target Gene'].str.split(' /// ', expand=True)[0]
        con = sqlite3.connect(fpath.replace('.txt', '.db'))
        df.to_sql('TargetScan', con, if_exists='replace', index=None)

    def TargetGene(self):
        fpath = os.path.join(self.root, 'database', 'feature_info.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'TargetScan'", con)
        df_grp = df.groupby('miRNA')

        contents = []
        for mir, df_sub in df_grp:
            tg = df_sub['Target Gene']
            contents.append([mir, ';'.join(tg)])
        df_res = pd.DataFrame(contents, columns=['miRNA', 'Target Gene'])
        df_res.to_excel('target_scan.xlsx', index=None)

    def mirTarbase(self):
        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'miRTarBase' WHERE miRNA LIKE '%hsa%'", con)
        df_grp = df.groupby('miRNA')

        contents = []
        for mir, df_sub in df_grp:
            tg = df_sub['Target Gene']
            sm = df_sub['Species_miRNA']
            ex = df_sub['Experiments']
            contents.append([mir.lower(), ';'.join(tg), ';'.join(ex)])
        df_res = pd.DataFrame(contents, columns=['miRNA', 'genes', 'Experiments'])
        df_res = df_res.sort_values(by=['miRNA'])
        df_res = df_res[df_res['genes'] > '']
        df_res.to_sql('miRTartBase_hsa', con, index=None, if_exists='replace')
        # df_res.to_excel('miRTartBase.xlsx', index=None)

    def target_genes(self):
        dirname = os.path.join(self.root, 'database', 'target_genes')

        out_con = sqlite3.connect(os.path.join(dirname, 'predictions_processed.db'))
        con = sqlite3.connect(os.path.join(dirname, 'predictions_processed_result.db'))

        for fname in ['miranda', 'rna22', 'ts']:
            fpath = os.path.join(dirname, 'predictions_processed_{}.txt'.format(fname))
            df = pd.read_csv(fpath, sep='\t', names=['miRNA', 'genes', 'start', 'end', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'])
            df = df[['miRNA', 'genes', 'start', 'end']]
            df.to_sql(fname, out_con, index=None, if_exists='replace')

            df_grp = df.groupby('miRNA')
            df_res = pd.DataFrame(index=df_grp.groups, columns=['genes'])
            df_res.index.name = 'miRNA'
            for mir, df_grp in df_grp:
                df_res.loc[mir, 'genes'] = ';'.join(df_grp['genes'])
            df_res.to_sql(fname, con, if_exists='replace')

    def get_intersection_all(self):
        dirname = os.path.join(self.root, 'database', 'target_genes')
        con = sqlite3.connect(os.path.join(dirname, 'predictions_processed_result.db'))

        dfs = {}
        mirs = []
        tlist = [x for x in Database.load_tableList(con) if '_pre' in x]
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='pre-miRNA')['genes']
            dfs[tname] = df
            mirs.append(set(df.index))

        mirs = set.intersection(*mirs)
        df_res = pd.DataFrame(index=mirs, columns=['genes'])
        for mir in mirs:
            genes = []
            for k, df in dfs.items():
                genes.append(set(df.loc[mir].split(';')))
            df_res.loc[mir, 'genes'] = ';'.join(sorted(list(set.intersection(*genes))))
        df_res.index.name = 'pre-miRNA'
        df_res.to_sql('common_all_pre', con, if_exists='replace')

    def get_intersection(self):
        dirname = os.path.join(self.root, 'database', 'target_genes')
        con = sqlite3.connect(os.path.join(dirname, 'predictions_processed_result.db'))

        dfs = {}
        mir = {}
        tlist = [x for x in Database.load_tableList(con) if '_pre' in x]
        for tname in tlist:
            dfs[tname] = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='pre-miRNA')
            mir[tname] = set(dfs[tname].index)

        source = {}
        mir_union = sorted(list(set.union(*mir.values())))
        for mi in mir_union:
            for tname, mlist in mir.items():
                if mi in mlist:
                    if mi not in source:
                        source[mi] = [tname]
                    else:
                        source[mi].append(tname)

        df_src = pd.DataFrame(index=source.keys(), columns=['source'])
        for k, v in source.items():
            df_src.loc[k, 'source'] = ','.join(v)

        df_res = pd.DataFrame(index=mir_union, columns=['genes'])
        df_res.loc[df_src.index, 'source'] = df_src['source']
        for mi in mir_union:
            genes = {}
            for tname, df in dfs.items():
                if mi not in df.index:
                    continue
                genes[tname] = list(set(df.loc[mi, 'genes'].split(';')))

            inter_genes = {}
            for tname, gene_list in genes.items():
                for gene in gene_list:
                    inter_genes[gene] = inter_genes.get(gene, 0) + 1

            genes = []
            for k, v in inter_genes.items():
                if v > 1:
                    genes.append(k)

            if len(genes) > 0:
                df_res.loc[mi, 'genes'] = ';'.join(genes)
        df_res = df_res.dropna(subset=['genes'])
        df_res.index.name = 'pre-miRNA'
        df_res.to_sql('common_pre', con, if_exists='replace')

    def get_union(self):
        dirname = os.path.join(self.root, 'database', 'target_genes')
        con = sqlite3.connect(os.path.join(dirname, 'predictions_processed_result.db'))

        dfs = {}
        mir = {}
        tlist = [x for x in Database.load_tableList(con) if '_pre' in x]
        for tname in tlist:
            dfs[tname] = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='pre-miRNA')
            mir[tname] = set(dfs[tname].index)

        source = {}
        mir_union = sorted(list(set.union(*mir.values())))
        for mi in mir_union:
            for tname, mlist in mir.items():
                if mi in mlist:
                    if mi not in source:
                        source[mi] = [tname]
                    else:
                        source[mi].append(tname)

        df_src = pd.DataFrame(index=source.keys(), columns=['source'])
        for k, v in source.items():
            df_src.loc[k, 'source'] = ','.join(v)

        df_res = pd.DataFrame(index=mir_union, columns=['genes'])
        df_res.loc[df_src.index, 'source'] = df_src['source']
        for mi in mir_union:
            genes = {}
            for tname, df in dfs.items():
                if mi not in df.index:
                    continue
                genes[tname] = list(set(df.loc[mi, 'genes'].split(';')))

            inter_genes = {}
            for tname, gene_list in genes.items():
                for gene in gene_list:
                    inter_genes[gene] = inter_genes.get(gene, 0) + 1

            genes = []
            for k, v in inter_genes.items():
                if v > 0:
                    genes.append(k)

            if len(genes) > 0:
                df_res.loc[mi, 'genes'] = ';'.join(genes)
        df_res = df_res.dropna(subset=['genes'])
        df_res.index.name = 'pre-miRNA'
        df_res.to_sql('union_pre', con, if_exists='replace')

    def targetscan(self):
        fpath = os.path.join(self.root, 'database', 'Predicted_Targets_Info.default_predictions.txt')
        df = pd.read_csv(fpath, sep='\t')
        df = df[df['Species ID'] == 9606].reset_index(drop=True)
        mir_fam = df['miR Family'].str.split('/')

        contents = []
        for idx in mir_fam.index:
            if idx % 1000 == 0 or idx + 1 == mir_fam.shape[0]:
                print('{:,d} / {:,d}'.format(idx, mir_fam.shape[0]))
            for i, val in enumerate(mir_fam.loc[idx]):
                if i > 0:
                    if val[0] != 'm':
                        val = 'mir-' + val
                contents.append([val.lower(), *df.loc[idx, ['Gene ID', 'Gene Symbol', 'Transcript ID']]])

        df = pd.DataFrame(contents, columns=['miRNA', 'Gene ID', 'Gene Symbol', 'Transcript ID'])
        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)

        df['miRNA'] = 'hsa-' + df['miRNA']
        df = df.sort_values(by=['miRNA'])
        df.to_sql('target_scan', con, index=None, if_exists='replace')

    def targetscan_grp(self):
        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)
        # df = pd.read_sql("SELECT * FROM 'target_scan'", con)
        # df_grp = df.groupby('miRNA')

        mir = pd.read_sql("SELECT DISTINCT miRNA FROM 'target_scan'", con)
        df_res = pd.DataFrame(index=mir['miRNA'], columns=['genes', 'transcript_id'])
        for i in mir.index:
            m = mir.loc[i, 'miRNA']
            df = pd.read_sql("SELECT * FROM 'target_scan' WHERE miRNA='{}'".format(m), con)
            df_res.loc[m, 'genes'] = ';'.join(df['Gene Symbol'])
            df_res.loc[m, 'transcript_id'] = ';'.join(df['Transcript ID'])
        df_res.to_sql('target_scan_grp', con, if_exists='replace')

    def targetscan_corr(self):
        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)

        fpath_tid = os.path.join(self.root, 'database/gencode', 'gene_tid.csv')
        df_tid = pd.read_csv(fpath_tid, sep='\t', index_col=1)

        fpath_corr = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100.db')
        con_corr = sqlite3.connect(fpath_corr)
        df_corr = pd.read_sql("SELECT * FROM 'corr'", con_corr, index_col='transcript_id')

        df = pd.read_sql("SELECT * FROM 'target_scan_grp'", con, index_col='miRNA')
        genes = df['genes'].str.split(';')
        mirlist = set.intersection(set(df.index), set(df_corr.columns))

        for mir in mirlist:
            corrs = []
            for gene in genes[mir]:
                if len(gene) == 0:
                    continue
                tids = df_tid.loc[gene, 'transcript_id']
                corr = df_corr.loc[set.intersection(set(df_corr.index), set(tids)), mir].max()
                if not math.isnan(corr):
                    corrs.append(corr)
            if not corrs:
                continue
            corrs = np.array(corrs)
            df.loc[mir, 'corr'] = ';'.join(corrs.round(4).astype(str))
            corrs_stats = np.array([corrs.mean(), corrs.std(), corrs.min(), corrs.max(), np.median(corrs)])
            df.loc[mir, 'corr_stats'] = ';'.join(corrs_stats.round(4).astype(str))

        df.to_sql('target_scan_grp', con, if_exists='replace')

    def get_amlan(self):
        fpath = os.path.join(self.root, 'database', 'important_genes_from_Amlan_each.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        out_con = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'predictions_processed_result.db'))
        contents = []
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='miRNA').astype(float)
            coeff_abs = df['coeff'].abs()
            coeff_abs = coeff_abs.sort_values(ascending=False)[:100]
            df = df.loc[coeff_abs.index]
            # coeff_mean = coeff_abs.mean()
            # coeff_std = coeff_abs.std()
            # df = df[coeff_abs > (coeff_mean + coeff_std)]
            genes = sorted(list(set(df.index)))
            contents.append([tname, ';'.join(genes), None, None])
        df = pd.DataFrame(contents, columns=['pre-miRNA', 'genes', 'coeff_mean', 'coeff_std'])
        df[['pre-miRNA', 'genes']].to_sql('amlan_pre', out_con, if_exists='replace', index=None)

    def test(self):
        fpath = os.path.join(self.root, 'database', 'important_genes_from_Amlan.xlsx')
        df_ref = pd.read_excel(fpath, index_col=0)
        df_ref = df_ref.dropna(subset=['genes'])

        df_high = pd.read_excel(os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_100.xlsx'))
        high_genes = set(df_high['gene_name'])

        for idx in df_ref.index:
            genes = set(df_ref.loc[idx, 'genes'].split(';'))
            compl = genes - high_genes
            print(len(compl))

    def comparison(self, hbw):
        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)

        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_{}_2.xlsx'.format(hbw))
        writer = pd.ExcelWriter(out_path, engine='xlsxwriter')
        for tname in ['miRTartBase_hsa', 'target_scan_grp']:
            df_ref = pd.read_sql("SELECT * FROM '{}' WHERE genes>''".format(tname), con, index_col='miRNA')

            df_high = pd.read_excel(os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.xlsx'.format(hbw)))
            high_genes = set(df_high['gene_name'])
            N = len(high_genes)

            fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}.db'.format(hbw))
            con_regr = sqlite3.connect(fpath)
            df = pd.read_sql("SELECT * FROM 'result'", con_regr, index_col='miRNA')
            mirs = set.intersection(set(df_ref.index), set(df.index))

            contents = []
            for mir in mirs:
                if isinstance(df.loc[mir, 'gene_name'], float):
                    continue

                genes = set(df.loc[mir, 'gene_name'].split(';'))
                if isinstance(df_ref.loc[mir, 'genes'], float):
                    continue

                genes_ref = set(df_ref.loc[mir, 'genes'].split(';'))
                genes_ref = [x for x in genes_ref if len(x) > 0]
                m = len(genes_ref)
                n = N - m
                q = len(set.intersection(genes, genes_ref))
                k = len(genes)
                contents.append([mir, m, n, q, k])
            pd.DataFrame(contents, columns=['mir', 'm', 'n', 'q', 'k']).to_excel(writer, sheet_name=tname, index=None)
        writer.save()
        writer.close()

    def phypher(self, hbw):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_{}_2.xlsx'.format(hbw))

        writer = pd.ExcelWriter(fpath, engine='xlsxwriter')
        for tname in ['miRTartBase_hsa', 'target_scan_grp']:
            df = pd.read_excel(fpath, index_col=0, sheet_name=tname)
            for idx in df.index:
                p = hypergeom.sf(df.loc[idx, 'q']-1, df.loc[idx, 'n'] + df.loc[idx, 'm'], df.loc[idx, 'm'], df.loc[idx, 'k'])
                df.loc[idx, 'phypher'] = p
                df.loc[idx, 'pvalue'] = np.log(1-p)
            df.to_excel(writer, sheet_name=tname)
        writer.save()
        writer.close()

    def plot(self, hbw):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'result_{}_2.xlsx'.format(hbw))
        tnames = ['miRTartBase_hsa', 'target_scan_grp']
        N = len(tnames)
        for i, tname in enumerate(tnames):
            df = pd.read_excel(fpath, index_col=0, sheet_name=tname)
            df['pvalue'] = df['pvalue'].replace(-np.inf, np.nan).replace(np.inf, np.nan)
            df = df.dropna(subset=['pvalue'])

            xaxis = range(df.shape[0])
            plt.subplot(N, 1, i+1)
            plt.scatter(xaxis, df['pvalue'])
            plt.ylabel('log(1-phypher)')
            plt.xticks(xaxis, df.index, fontsize=6, rotation=30)
            plt.title('max: {:0.2f}, min: {:0.2f}, mean: {:0.2f}'.format(df['pvalue'].max(), df['pvalue'].min(), df['pvalue'].mean()))
        plt.show()
        # plt.savefig(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'log_phypher.png'))

    def matrue_to_pre(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'mirna_primirna_mapping.txt')
        df = pd.read_csv(fpath, sep='\t')

        con = sqlite3.connect(fpath.replace('.txt', '.db'))
        df.to_sql('mirna_primirna_mapping', con, if_exists='replace', index=None)


if __name__ == '__main__':
    mg = mir_gene()
    # mg.test()
    if mg.hostname == 'mingyu-Precision-Tower-781':
        # mg.targetscan()
        # mg.mirTarbase()
        # mg.targetscan_grp()
        # mg.targetscan_corr()
        mg.to_server()
    else:
        # mg.matrue_to_pre()
        # mg.get_intersection_all()
        # mg.get_intersection()
        mg.get_union()
        # mg.target_genes()
        # mg.get_amlan()
        # mg.mirTarbase()
        # mg.targetscan()
        # mg.targetscan_grp()

        # mg.comparison(100)
        # mg.phypher(100)
        # mg.plot(100)
