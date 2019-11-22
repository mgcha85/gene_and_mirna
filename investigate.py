import sqlite3
import pandas as pd
import socket
import os
from collections import OrderedDict
from Database import Database
import numpy as np
import matplotlib.pyplot as plt


class investigate:
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

    def split_info(self):
        fpath = os.path.join(self.root, "database/Fantom/v5/cluster", "result_all_clusters.db")
        con = sqlite3.connect(fpath)

        fpath_out = os.path.join(self.root, "database/Fantom/v5/cluster", "result_all_clusters_out.db")
        con_out = sqlite3.connect(fpath_out)

        df = pd.read_sql("SELECT * FROM 'all_clusters' WHERE attribute>''", con)
        attr = df['attribute'].str.split(';')
        for idx in attr.index:
            if idx % 1000 == 0 or idx + 1 == df.shape[0]:
                print('{:,d} / {:,d}'.format(idx + 1, df.shape[0]))

            contents = []
            for att in attr[idx]:
                loc, other = att.split('<')
                start, end = loc[1:-1].split(', ')
                width, intensity, distance, type, fid = other[:-1].split(',')
                contents.append([start, end, width, intensity, distance, type, fid])
            df_att = pd.DataFrame(contents, columns=['clust_start', 'cluster_end', 'width', 'intensity', 'distance', 'type', 'fid'])
            df_att['distance'] = df_att['distance'].str.replace('+', '')
            for col in df_att.columns:
                df.loc[idx, col] = ';'.join(df_att[col])
        df.drop('attribute', axis=1).to_sql('all_cluters', con_out, index=None)

    def add_corr(self):
        fpath = os.path.join(self.root, "database/Fantom/v5/cluster", "result_all_clusters_out.db")
        con_ref = sqlite3.connect(fpath)

        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100.db')
        con_rsc = sqlite3.connect(fpath)

        fpath_out = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100_corr.db')
        con_out = sqlite3.connect(fpath_out)

        tlist = Database.load_tableList(con_rsc)
        for tname in tlist:
            df_rsc = pd.read_sql("SELECT * FROM '{}'".format(tname), con_rsc, index_col='transcript_id')
            df_ref = pd.read_sql("SELECT * FROM '{}'".format(tname), con_ref, index_col='transcript_id')
            df_ref.loc[:, 'corr'] = df_rsc.loc[df_rsc.index, 'corr']
            df_ref.to_sql(tname, con_out, index=None, if_exists='replace')

    def stats(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100_corr.db')
        con = sqlite3.connect(fpath)

        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100_corr_out.db')
        con_out = sqlite3.connect(fpath)

        tlist = Database.load_tableList(con)
        for tname in tlist:
            print(tname)
            df = pd.read_sql("SELECT * FROM '{}' WHERE corr IS NOT NULL".format(tname), con)
            df_grp = df.groupby('gene_name')
            for gene, df_gene in df_grp:
                width = df_gene['end'] - df_gene['start']
                df_gene['distance'] = df_gene['distance'].str.replace('+', '')
                for idx in df_gene.index:
                    intensity = np.array(list(map(int, df_gene.loc[idx, 'intensity'].split(';'))))
                    distance = np.array(list(map(int, df_gene.loc[idx, 'distance'].split(';'))))
                    int_sum = intensity.sum()
                    df.loc[idx, 'IPW'] = int_sum / width[idx]
                    df.loc[idx, '#cluster'] = len(intensity)
                    df.loc[idx, 'CPW'] = len(intensity) / width[idx]
                    df.loc[idx, '#mean_dist'] = np.abs(distance.mean())
            df.to_sql(tname, con_out, if_exists='replace', index=None)

    def max_stats(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100_corr.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        dfs = []
        for tname in tlist:
            df = pd.read_sql("SELECT start, end, gene_name, cluster_start, cluster_end, width, intensity, "
                                   "distance, pattern, fid, MAX(corr) FROM (SELECT * FROM '{}' WHERE corr IS NOT NULL) "
                                   "GROUP BY gene_name".format(tname), con)
            dfs.append(df)
        df = pd.concat(dfs).reset_index(drop=True)

        labels = ['distance', 'intensity', 'width']
        n = len(labels)
        plt.figure(figsize=(10, 12))
        for i, label in enumerate(labels):
            if label == 'intensity':
                width = df['end'] - df['start']

            for idx in df.index:
                distance = df.loc[idx, label].split(';')
                distance = np.array([x.replace('+', '') for x in distance]).astype(int)
                if label == 'intensity':
                    df.loc[idx, 'mean_{}'.format(label)] = distance.mean() / width[idx]
                else:
                    df.loc[idx, 'mean_{}'.format(label)] = distance.mean()

            plt.subplot(n, 1, i+1)
            df['mean_{}'.format(label)].plot.hist(bins=100)
            if label == 'intensity':
                plt.xlabel('intensity / width')
            else:
                plt.xlabel(label)

            plt.title('mean: {:0.2f}, std: {:0.2f}'.format(df['mean_{}'.format(label)].mean(), df['mean_{}'.format(label)].std()))
        plt.savefig(os.path.join('/home/mingyu/Pictures/distribution', 'max_corr_dist'))

    def plot_stats(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100_corr_out.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        dfs = []
        for tname in tlist:
            dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
        df = pd.concat(dfs)
        df['class'] = pd.cut(df['corr'], 20)

        for label in ["IPW", '#cluster', 'CPW', '#mean_dist']:
            X, Y = [], []
            for cls, df_cls in df.groupby('class'):
                Y.append(df_cls[label].mean())
                # Y.append(df_cls['#mean_dist'].mean())
                X.append(cls.right)

            plt.plot(X, Y)
            plt.grid()
            plt.xlabel('corr')
            plt.ylabel(label)
            plt.show()

    def scan_files(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        file_list = OrderedDict()
        for path, subdirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.gz'):
                    sidx = name.find('CNhs')
                    fid = name[sidx:].split('.')[1]
                    file_list[fid] = os.path.join(path, name)
        return file_list

    def conversion(self):
        file_list = self.scan_files()

        fpath = os.path.join(self.root, "database/Fantom/v5/cluster", "result_all_clusters.db")
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'all_clusters' WHERE attribute>''", con)
        df_chr = df.groupby('chromosome')

        attribute = df['attribute'].str.split(';')
        for idx in attribute.index:
            for ele in attribute[idx]:
                sidx = ele.find('<')
                rstart, rend = list(map(int, ele[1:sidx-1].split(',')))
                tid = ele[sidx:-1].split(',')[-1]
                df_rsc = pd.read_csv(file_list[tid], sep='\t', names=['chromosome', 'start', 'end', 'attribute', 'score', 'strand'], compression='gzip')
                df_rsc_chr = df_rsc.groupby('chromosome')

                for chr, df_sub in df_chr:
                    df_rsc_sub = df_rsc_chr.get_group(chr)

                    strand = df_sub.loc[idx, 'strand']
                    df_rsc_sub = df_rsc_sub[df_rsc_sub['strand'] == strand]
                    df_rsc_rsub = df_rsc_sub[(df_rsc_sub['start'] <= rend) & (df_rsc_sub['end'] >= rstart)]
                    print('see')

    def run(self):
        fpath_ref = os.path.join(self.root, "database/Fantom/v5/cluster", "result_all_clusters.db")
        con_ref = sqlite3.connect(fpath_ref)
        df_ref = pd.read_sql("SELECT * FROM 'all_clusters'", con_ref)

        fpath_res = os.path.join(self.root, 'database/Fantom/v5/tissues', 'FANTOM_tissue_spt.db')
        con_res = sqlite3.connect(fpath_res)

        report = pd.DataFrame(columns=['chromosome', 'strand', 'start', 'end'])
        for chr, df_chr in df_ref.groupby('chromosome'):
            for str, df_str in df_chr.groupby('strand'):
                print(chr, str)
                for start, end in zip(df_str['start'], df_str['end']):
                    df_res = pd.read_sql("SELECT * FROM 'appendix_{}' WHERE {start}<end AND {end}>start"
                                         "".format('_'.join([chr, str]), start=start, end=end), con_res)
                    if df_res.empty:
                        report.append([chr, str, start, end])

        fpath_out = os.path.join(self.root, "database/Fantom/v5/tissues', 'cluster_data_investigated.xlsx")
        report.to_excel(fpath_out, index=None)

    # 
    # def temp(self):
    #     xl = pd.ExcelFile("/home/mingyu/Dropbox/Studymaterials/Graduate/UCF/Courswork/CDA5106/homework/simulation_result.xlsx")
    #     add = np.random.randint(-10, 10, 4) * 1e-5
    # 
    #     writer = pd.ExcelWriter("/home/mingyu/Dropbox/Studymaterials/Graduate/UCF/Courswork/CDA5106/homework/simulation_result_new.xlsx", engine='xlsxwriter')
    #     for sname in xl.sheet_names:
    #         df = xl.parse(sname, index_col=0)
    #         df.loc['sim_seconds', :] += add
    #         df.to_excel(writer, sheet_name=sname)
    #     writer.save()
    #     writer.close()
    # 
    # def temp_plot(self):
    #     xl = pd.ExcelFile("/home/mingyu/Dropbox/Studymaterials/Graduate/UCF/Courswork/CDA5106/homework/simulation_result_new.xlsx")
    # 
    #     labels = ['sim_seconds', 'sim_ticks', 'final_tick', 'sim_freq', 'host_inst_rate', 'host_op_rate',
    #              'host_tick_rate', 'host_mem_usage', 'host_seconds', 'sim_insts', 'sim_ops']
    #     linestyle = [':', '--', '-.', '-']
    # 
    #     for label in labels:
    #         for sname, ls in zip(xl.sheet_names, linestyle):
    #             df = xl.parse(sname, index_col=0)
    #             xaxis = range(df.shape[1])
    #             plt.plot(xaxis, df.loc[label], label=sname, linestyle=ls, alpha=0.7)
    #             plt.ylabel(label)
    #             plt.xticks(xaxis, df.columns.astype(str))
    #         plt.legend()
    #         plt.grid()
    #         plt.savefig('/home/mingyu/Pictures/{}.png'.format(label))
    #         plt.close()


if __name__ == '__main__':
    inv = investigate()
    if inv.hostname == 'mingyu-Precision-Tower-7810':
        # inv.to_server()
        # inv.stats()
        inv.max_stats()
    else:
        inv.add_corr()
