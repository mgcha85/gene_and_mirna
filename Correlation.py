import socket
import os
import sqlite3
import pandas as pd
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import pickle as pkl
from Database import Database
from Server import Server
import sys



class Correlation:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
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

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def add_type(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df_fantom = pd.read_csv(fpath)

        fpath_fan = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath_fan)

        contents = []
        for idx in df_fantom.index:
            if idx % 1000 == 0 or idx == df_fantom.shape[0] - 1:
                print('{:0.2f}%'.format(100 * (idx + 1) / df_fantom.shape[0]))

            chromosome = df_fantom.loc[idx, 'chromosome']
            strand = df_fantom.loc[idx, 'strand']
            start = df_fantom.loc[idx, 'start']
            end = df_fantom.loc[idx, 'end']

            df_mir = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE chromosome='{}' AND "
                                       "strand='{}' AND tss>={start} AND tss<{end}"
                                       "".format(chromosome, strand, start=start, end=end), con)

            df_gene = pd.read_sql_query("SELECT * FROM 'human_gene' WHERE chromosome='{}' AND start>={start} AND "
                                        "start<{end}".format(chromosome, strand, start=start, end=end), con)
            if not df_mir.empty:
                mirna = ';'.join(df_mir['premiRNA'])
                contents.append(['miRNA', mirna])
            elif not df_gene.empty:
                gene = ';'.join(df_gene['gene'])
                contents.append(['gene', gene])
            else:
                contents.append(['other', None])

        df_type = pd.DataFrame(contents, columns=['tss-type', 'name'])
        df_fantom = pd.concat([df_fantom, df_type], axis=1)
        df_fantom.to_csv(fpath.replace('.csv', '2.csv'), index=None)

    def load_reference(self, tissue):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.final_confirmed_tss.db')
        con = sqlite3.connect(fpath)
        if Database.checkTableExists(con, tissue):
            return pd.read_sql_query("SELECT * FROM '{}'".format(tissue), con)

    def run(self):
        RANGE = 500
        out_columns = ['chromosome', 'start', 'end', 'strand', 'FPKM (fantom)', 'TPM (fantom)', 'FPKM (rna)',
                       'TPM (rna)', 'replication']

        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue.db')
        con_fan = sqlite3.connect(fpath)

        fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_tissue.db')
        con_rna = sqlite3.connect(fpath)

        tlist_fan = Database.load_tableList(con_fan)
        tlist_rna = Database.load_tableList(con_rna)

        con_out = sqlite3.connect(os.path.join(self.root, 'Papers/complement', 'correlation.db'))

        reports = {}
        for tname_fan in tlist_fan:
            print(tname_fan)
            df_ref = self.load_reference(tname_fan)
            if df_ref is None:
                continue

            # df_ref = df_ref[:1000]
            contents = []
            for idx in df_ref.index:
                if idx % 100 == 0 or idx + 1 == df_ref.shape[0]:
                    print('{:0.2f}%'.format(100 * (idx + 1) / df_ref.shape[0]))
                strand = df_ref.loc[idx, 'strand']
                if strand == '+':
                    tss = df_ref.loc[idx, 'start']
                else:
                    tss = df_ref.loc[idx, 'end']

                start = tss - RANGE
                end = tss + RANGE
                chromosome = df_ref.loc[idx, 'chromosome']

                df_fan = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND NOT "
                                           "start>{end} AND NOT end<{start}".format(tname_fan, chromosome, start=start,
                                                                                    end=end), con_fan)

                tnames_rna = [x for x in tlist_rna if tname_fan in x]
                for tname_rna in tnames_rna:
                    df_rna = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND NOT "
                                               "start>{end} AND NOT end<{start}".format(tname_rna, chromosome,
                                                                                        start=start, end=end), con_rna)
                    contents.append([chromosome, start, end, strand, df_fan['FPKM'].astype(float).sum(), df_fan['TPM'].astype(float).sum(),
                                     df_rna['FPKM'].astype(float).sum(), df_rna['TPM'].astype(float).sum(),
                                     tname_rna])
            reports[tname_fan] = pd.DataFrame(data=contents, columns=out_columns)

        for tissue, df in reports.items():
            df.to_sql(tissue, con_out, index=None, if_exists='replace')

    def split(self):
        fname = 'correlation_report'
        df = pd.read_excel('{}.xlsx'.format(fname), index_col=0)
        df.index.name = 'miRNA'

        N = 500  # number of split tables

        out_path = os.path.join(self.root, 'database', '{}.db'.format(fname))
        con = sqlite3.connect(out_path)

        M = int((df.shape[1] + N - 1) / N)
        k = len(str(M))
        for i in range(M):
            idx = i * N
            df_spt = df.iloc[:, idx: idx + N]
            df_spt.to_sql(fname + '_{}'.format(str(i).zfill(k)), con, if_exists='replace')

    def figure(self):
        df = pd.read_excel('correlation_report.xlsx', index_col=0)

        fig = plt.figure()
        # ax = Axes3D(fig)
        ax1 = fig.add_subplot(111, projection='3d')

        x = range(df.shape[1])
        y = range(df.shape[0])
        X, Y = np.meshgrid(x, y)
        values = df.values.flatten()
        Z = np.zeros_like(values)

        # color_values run= plt.cm.jet(df.values.tolist())
        # ax1.bar3d(X, Y, df.values, dx=1, dy=1, dz=1, color=color_values)
        ax1.bar3d(X.flatten(), Y.flatten(), Z, dx=1, dy=1, dz=values)
        plt.savefig('correlation_report.png')

    def filter(self):
        from Database import Database

        fpath = os.path.join(self.root, 'database', 'correlation_report.db')
        con = sqlite3.connect(fpath)
        tableList = Database.load_tableList(con)
        thres = 0.8

        dfs = {'GENE': [], 'miRNA': [], 'value': []}
        for tname in tableList:
            print(tname)
            df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con, index_col='miRNA')
            df = df[~df.index.duplicated(keep='first')]

            for idx in df.index:
                for col in df.columns:
                    val = df.loc[idx, col]
                    # if not isinstance(val, float):
                    #     print('see')
                    if abs(val) > thres:
                        dfs['miRNA'].append(idx)
                        dfs['GENE'].append(col)
                        dfs['value'].append(val)
        df_res = pd.DataFrame(dfs)
        df_res.to_excel('correlation_report2_{:.1f}.xlsx'.format(thres), index=None)

    def profile_gene_mirna(self):
        fname = 'correlation_report2_0.9.xlsx'
        df = pd.read_excel(fname, index_col=0)

        con = sqlite3.connect(os.path.join(self.root, 'database', 'fantom5.db'))
        df_gene = pd.read_sql_query("SELECT * FROM 'human_gene'", con, index_col='gene')

        report = []
        for i, gene in enumerate(df.index):
            if i % 100 == 0 or i + 1 == df.shape[0]:
                print('{:.2f}%'.format(100 * (i + 1) / df.shape[0]))

            if ';' in gene:
                continue
            if gene not in df_gene.index:
                continue

            gene_chromosome = df_gene.loc[gene, 'chromosome']
            gene_start = df_gene.loc[gene, 'start']
            gene_end = df_gene.loc[gene, 'end']
            mirna = df.loc[gene, 'miRNA']

            for mir in mirna:
                df_mir = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE premiRNA='{}'".format(mir), con)
                mir_gene = ';'.join(df_mir['Gene'].astype(str))
                mirna_chromosome = ';'.join(df_mir['chromosome'])
                mirna_start = ';'.join(df_mir['start'].astype(str))
                mirna_end = ';'.join(df_mir['end'].astype(str))
                mirna_tss = ';'.join(df_mir['tss'].astype(str))
                mirna_strand = ';'.join(df_mir['strand'])

                report.append([gene, gene_chromosome, gene_start, gene_end, mir_gene, mir, mirna_chromosome, mirna_start, mirna_end, mirna_tss, mirna_strand])

        df_rep = pd.DataFrame(data=report, columns=['gene', 'gene_chromosome', 'gene_start', 'gene_end', 'mir_gene', 'mir', 'mirna_chromosome', 'mirna_start', 'mirna_end', 'mirna_tss', 'mirna_strand'])
        df_rep.to_excel(fname, index=None)

    def loc_info_analysis(self):
        fname = 'correlation_loc_info_0.9.xlsx'
        df = pd.read_excel(fname)
        for idx in df.index:
            gene_chromosome = df.loc[idx, 'gene_chromosome']
            mirna_chromosome = df.loc[idx, 'mirna_chromosome']

            gene = df.loc[idx, 'gene']
            mir_gene = df.loc[idx, 'mir_gene']

            mirna_tss = df.loc[idx, 'mirna_tss']
            mirna_strand = df.loc[idx, 'mirna_strand']

            if gene == mir_gene:
                df.loc[idx, 'same gene'] = 'yes'
                df.loc[idx, 'same chrom'] = 'yes'

                if mirna_strand == '+':
                    df.loc[idx, 'distance'] = abs(mirna_tss - df.loc[idx, 'gene_start'])
                else:
                    df.loc[idx, 'distance'] = abs(mirna_tss - df.loc[idx, 'gene_end'])

            elif gene_chromosome == mirna_chromosome:
                df.loc[idx, 'same chrom'] = 'yes'
                df.loc[idx, 'same gene'] = 'no'
                if mirna_strand == '+':
                    df.loc[idx, 'distance'] = abs(mirna_tss - df.loc[idx, 'gene_start'])
                else:
                    df.loc[idx, 'distance'] = abs(mirna_tss - df.loc[idx, 'gene_end'])
            else:
                df.loc[idx, 'same gene'] = 'no'
                df.loc[idx, 'same chrom'] = 'no'
        
        df.to_excel(fname, index=None)

    def stats(self):
        df = pd.read_excel('correlation_loc_info_0.9.xlsx')
        df_grp = df.groupby('same chrom')

        report = {}
        for tf, df_sub in df_grp:
            report[tf] = df_sub.shape[0]

        df_sub = df_grp.get_group('yes')
        ele_hist, ele_bin = np.histogram(df_sub['distance'], bins=10)
        print(ele_hist)
        xaxis = range(len(ele_hist))
        plt.bar(xaxis, ele_hist)
        plt.xticks(xaxis, ele_bin[1:], rotation=30, fontsize=6)
        plt.show()


if __name__ == '__main__':
    cor = Correlation()
    if cor.hostname == 'mingyu-Precision-Tower-781':
        cor.to_server()
        # cor.run()
    else:
        cor.run()
