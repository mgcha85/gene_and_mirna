import pandas as pd
import sqlite3
import socket
import os


class Report:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def get_sample_fan(self):
        columns = ['chromosome', 'start', 'end', 'loc', 'score', 'strand']

        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues/adipose')
        flist = os.listdir(dirname)

        dfs = []
        idx = []
        for fname in flist:
            fpath = os.path.join(dirname, fname)
            for df_chunk in pd.read_csv(fpath, sep='\t', names=columns, compression='gzip', chunksize=1e3, index_col=3):
                dfs.append(df_chunk)
                idx.append(set(dfs[-1].index))
                break

        idx = sorted(list(set.union(*idx)))
        df_res = pd.DataFrame(index=idx, columns=['chromosome', 'start', 'end', 'avg_score', 'strand'])

        cols = []
        for i, df in enumerate(dfs):
            cols.append('score (Rep{})'.format(i))
            df_res.loc[df.index, ['chromosome', 'start', 'end', 'strand']] = df[['chromosome', 'start', 'end', 'strand']]
            df_res.loc[df.index, cols[-1]] = df['score']

        df_res['avg_score'] = df_res[cols].mean(axis=1)
        df_res.to_excel(os.path.join(dirname, 'avg_score.xlsx'), index=None)

    def get_sample_rna(self):
        df = pd.read_excel('RNA-seq_data_structure.xlsx')
        columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'gene_name',
                   'transcript_id']

        dirname = os.path.join(self.root, 'database/RNA-seq/out')
        fpath = os.path.join(dirname, 'RNA_seq.db')
        con = sqlite3.connect(fpath)
        for src, df_src in df.groupby('tissue'):
            if src != 'fat':
                continue

            dfs = []
            index = []
            for fid in df_src['fid']:
                df_fid = pd.read_sql_query("SELECT * FROM '{}' WHERE FPKM>''".format(fid), con)[:1000]
                df_fid.to_csv(os.path.join(dirname, fid + '.csv'), index=None)
                df_fid.index = df_fid['chromosome'] + ':' + df_fid['start'].astype(str) + '-' + df_fid[
                    'end'].astype(str)
                index.append(set(df_fid.index))
                dfs.append(df_fid)

            index = sorted(list(set.union(*index)))
            df_res = pd.DataFrame(index=index, columns=df_fid.columns)
            df_res[['cov', 'FPKM', 'TPM']] = 0
            for i, df_fid in enumerate(dfs):
                df_res.loc[df_fid.index, ['cov', 'FPKM', 'TPM']] += df_fid[['cov', 'FPKM', 'TPM']].astype(float)
                df_res.loc[df_fid.index, 'FPKM (Rep {})'.format(i + 1)] = df_fid['FPKM'].astype(float)
                df_res.loc[df_fid.index, columns] = df_fid[columns]

            df_res[['cov', 'FPKM', 'TPM']] /= len(dfs)
            df_res.to_excel(os.path.join(dirname, 'avg_score.xlsx'), index=None)

    def scores_by_tissues(self):
        from Database import Database

        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna.db')
        con = sqlite3.connect(fpath)

        writers = {'fantom': pd.ExcelWriter(fpath.replace('.db', '_fan.xlsx'), engine='xlsxwriter'),
                   'rna-seq': pd.ExcelWriter(fpath.replace('.db', '_rna.xlsx'), engine='xlsxwriter')}
        tlist = Database.load_tableList(con)
        for tname in tlist:
            source, chr, str = tname.split('_')
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            df.to_excel(writers[source], sheet_name=tname, index=None)

        for source, writer in writers.items():
            writer.save()
            writer.close()

    def to_excel(self):
        from Database import Database
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out')
        fpath = os.path.join(dirname, 'regression_100.db')
        con = sqlite3.connect(fpath)

        tlist = Database.load_tableList(con)
        writer = pd.ExcelWriter(os.path.join(fpath.replace('.db', '.xlsx')), engine='xlsxwriter')
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            df.to_excel(writer, sheet_name=tname, index=None)
        writer.save()
        writer.close()


if __name__ == '__main__':
    rep = Report()
    rep.to_excel()
    # rep.scores_by_tissues()
