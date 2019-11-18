import sqlite3
import pandas as pd
import socket
import os
from collections import OrderedDict


class investigate:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def split_info(self):
        fpath = os.path.join(self.root, "database/Fantom/v5/cluster", "result_all_clusters.db")
        con = sqlite3.connect(fpath)

        df = pd.read_sql("SELECT * FROM 'all_clusters' WHERE attribute>''", con)
        attr = df['attribute'].str.split(';')
        contents = []
        for idx in attr.index:
            if idx % 1000 == 0 or idx + 1 == df.shape[0]:
                print('{:,d} / {:,d}'.format(idx + 1, df.shape[0]))
            for att in attr[idx]:
                loc, other = att.split('<')
                start, end = loc[1:-1].split(', ')
                width, intensity, distance, type, fid = other[:-1].split(',')
                contents.append([start, end, width, intensity, distance, type, fid])
            df_att = pd.DataFrame(contents, columns=['clust_start', 'cluster_end', 'width', 'intensity', 'distance', 'type', 'fid'])
            df_att['distance'] = df_att['distance'].str.replace('+', '')
            for col in df_att.columns:
                df.loc[idx, col] = ';'.join(df_att[col])

        fpath = os.path.join(self.root, "database/Fantom/v5/cluster", "result_all_clusters_out.db")
        df.drop('attribute', axis=1).to_sql('all_cluters', sqlite3.connect(fpath), index=None)

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


if __name__ == '__main__':
    inv = investigate()
    inv.split_info()
