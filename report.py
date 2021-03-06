import pandas as pd
import sqlite3
import socket
import os


class Report:
    def __init__(self, root):
        self.root = root

    def get_lasso_result(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out')
        con = sqlite3.connect(os.path.join(dirname, 'regression_100_nz.db'))
        df = pd.read_sql("SELECT miRNA, gene_name FROM 'result'".format(), con, index_col='miRNA')
        gene_name = df['gene_name'].str.split(';')

        n = []
        for gn in gene_name:
            n.append(len(gn))
        df['n_genes'] = n
        df.to_excel(os.path.join(dirname, 'regression_100_nz.xlsx'))

    def compare_others(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out')
        fpath = os.path.join(dirname, 'regression_100_nz.xlsx')

        df = pd.read_excel(fpath, index_col=0)

        con_ref = sqlite3.connect(os.path.join(self.root, 'database', 'target_genes.db'))
        df_ts = pd.read_sql("SELECT * FROM 'target_scan_grp'", con_ref)
        df_mt = pd.read_sql("SELECT * FROM 'miRTartBase_hsa'", con_ref)

        df_res = pd.DataFrame(index=df.index)
        for mir in df.index:
            lgenes = df.loc[mir, 'gene_name'].split(';')
            tidx = df_ts[df_ts['pre-miRNA'] == mir].index

            intersection = None
            max = 0
            for t in tidx:
                genes = df_ts.loc[t, 'genes'].split(';')
                inter = set.intersection(set(lgenes), set(genes))
                if max < len(inter):
                    intersection = inter
            if intersection:
                df_res.loc[mir, '∩ ts'] = ';'.join(intersection)
                df_res.loc[mir, '# ∩ ts'] = len(intersection)

            tidx = df_mt[df_mt['pre-miRNA'] == mir].index

            intersection = None
            max = 0
            for t in tidx:
                genes = df_mt.loc[t, 'genes'].split(';')
                inter = set.intersection(set(lgenes), set(genes))
                if max < len(inter):
                    intersection = inter
            if intersection:
                df_res.loc[mir, '∩ mt'] = ';'.join(intersection)
                df_res.loc[mir, '# ∩ mt'] = len(intersection)
        df_res.to_excel(os.path.join(dirname, 'regression_100_nz2.xlsx'))

    def get_compare_cell_tissue(self):
        dfs = {}
        for type in ['tissues', 'cell_lines']:
            fpath = os.path.join(self.root, 'database/Fantom/v5/{}/out'.format(type), 'regression_100_nz.db')
            con = sqlite3.connect(fpath)
            dfs[type] = pd.read_sql("SELECT miRNA, gene_name FROM 'result'", con, index_col='miRNA')

        df_res = pd.DataFrame(index=dfs['cell_lines'].index)
        for mir in dfs['cell_lines'].index:
            genes = []
            for type in ['tissues', 'cell_lines']:
                genes.append(set(dfs[type].loc[mir, 'gene_name'].split(';')))

            intersection = set.intersection(*genes)
            union = set.union(*genes)
            df_res.loc[mir, '#tissue'] = len(genes[0])
            df_res.loc[mir, '#cell_lines'] = len(genes[1])
            df_res.loc[mir, '#∩'] = len(intersection)
            df_res.loc[mir, '#U'] = len(union)

        df_res['#∩/#U'] = df_res['#∩'] / df_res['#U']
        df_res['#∩/#small'] = df_res['#∩'] / df_res[['#tissue', '#cell_lines']].min(axis=1)
        df_res['#∩/#large'] = df_res['#∩'] / df_res[['#tissue', '#cell_lines']].max(axis=1)
        df_res['larger'] = df_res[['#∩/#small', '#∩/#large']].max(axis=1)
        df_res.to_excel(os.path.join(self.root, 'database/Fantom/v5/{}/out'.format(type), 'get_compare_cell_tissue.xlsx'))

    def missing_mir(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz.db')
        con = sqlite3.connect(fpath)
        df_corr = pd.read_sql("SELECT * FROM 'corr'", con, index_col='transcript_id')
        df_x = pd.read_sql("SELECT * FROM 'X'", con, index_col='tid')
        df_y = pd.read_sql("SELECT * FROM 'Y'", con, index_col='miRNA')
        df_res = pd.read_sql("SELECT * FROM 'result'", con, index_col='miRNA')

        missing_mir = set(df_corr.columns) - set(df_res.index)
        df_report = pd.Series(index=missing_mir)
        for m in missing_mir:
            df_report[m] = df_y.loc[m].sum()
        df_report.to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'missing_mir.xlsx'))

    def stats_mannwhitneyu(self):
        fpath = os.path.join(self.root, "database/Fantom/v5/cell_lines/out", 'mannwhitneyu.xlsx')
        df = pd.read_excel(fpath, index_col=0)
        df_stats = pd.DataFrame(index=['ts', 'mi'])

        for lab in df_stats.index:
            df_sub = df.dropna(subset=['stat ({})'.format(lab)])
            df_stats.loc[lab, '# <0.01'] = df_sub[df_sub['p-value ({})'.format(lab)] < 0.01].shape[0]
            df_stats.loc[lab, '# <0.001'] = df_sub[df_sub['p-value ({})'.format(lab)] < 0.001].shape[0]

            df_sub2 = df_sub[df_sub['# genes ({})'.format(lab)] >= 10]
            df_stats.loc[lab, '# <0.01 (genes>=10)'] = df_sub2[df_sub2['p-value ({})'.format(lab)] < 0.01].shape[0]
            df_stats.loc[lab, '# <0.001 (genes>=10)'] = df_sub2[df_sub2['p-value ({})'.format(lab)] < 0.001].shape[0]
        print(df_stats)

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
                df_fid = pd.read_sql("SELECT * FROM '{}' WHERE FPKM>''".format(fid), con)[:1000]
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

    def get_sample_corr(self):
        from Database import Database
        import numpy as np

        np.random.seed(0)
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        dfs = []
        for tname in tlist:
            chromosome, strand = tname.split('_')
            df = pd.read_sql("SELECT * FROM '{}' WHERE transcript_type='protein_coding' AND corr>0.9".format(tname), con)
            if df.empty:
                continue
            df.loc[:, 'chromosome'] = chromosome
            df.loc[:, 'strand'] = strand
            dfs.append(df)
        df = pd.concat(dfs).reset_index(drop=True)

        ridx = np.random.randint(0, df.shape[0], 3)
        df = df.iloc[ridx]

        fpath_score = os.path.join(self.root, 'database/Fantom/v5/tissues', 'sum_fan_rna_100.db')
        con_score = sqlite3.connect(fpath_score)

        contents = []
        for idx in df.index:
            chromosome = df.loc[idx, 'chromosome']
            strand = df.loc[idx, 'strand']
            start = df.loc[idx, 'start']
            end = df.loc[idx, 'end']
            df_fan = pd.read_sql("SELECT * FROM 'fantom_{}_{}' WHERE start={} AND end={}".format(chromosome, strand, start, end), con_score)
            df_rna = pd.read_sql("SELECT * FROM 'rna-seq_{}_{}' WHERE start={} AND end={}".format(chromosome, strand, start, end), con_score)
            contents.append([chromosome, strand, start, end, ','.join(df_fan.iloc[0, 2:].round(4).astype(str)), ','.join(df_rna.iloc[0, 2:].round(4).astype(str))])
        pd.DataFrame(data=contents, columns=['chromosome', 'strand', 'start', 'end', 'fantom', 'rna-seq']).to_excel('3_samples.xlsx', index=None)

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

    def to_sql_each(self):
        fpath = '/home/mingyu/Bioinformatics/database/important_genes_from_Amlan.txt'
        df = pd.read_csv(fpath, sep='\t', names=['miRNA', 'genes', 'coeff'])
        genes = df['genes'].str.split(';')
        corrs = df['coeff'].str.split(';')

        con = sqlite3.connect(fpath.replace('.txt', '.db'))
        for idx in genes.index:
            print(idx + 1, genes.shape[0])
            if isinstance(genes[idx], float):
                continue

            gs = []
            for i, col in enumerate(genes[idx]):
                if ' /// ' in col:
                    gs.append(col.split(' /// ')[0])
                else:
                    gs.append(col)

            ser = pd.DataFrame(data=[corrs[idx]], columns=gs, index=['coeff']).T
            ser.index.name = 'gene'
            ser.to_sql(df.loc[idx, 'miRNA'], con, if_exists='replace')

    # def to_sql(self):
    #     fpath = '/home/mingyu/Bioinformatics/database/human_cell_line_hCAGE_100_score_corr2_with_target_indicator_ts.csv'
    #     df = pd.read_csv(fpath)
    #     df = df.rename(columns={'GENEs': 'genes_org'})
    #     ti = df['Target Indicator'].str.split(';')
    #     genes = df['genes_org'].str.split(';')
    #
    #     for idx in ti.index:
    #         tr, gr = ti[idx], genes[idx]
    #
    #         gs = []
    #         for t, g in zip(tr, gr):
    #             if t == '1':
    #                 gs.append(g)
    #
    #         gs = ';'.join(gs)
    #         df.loc[idx, 'genes'] = gs
    #
    #     df = df[df['genes'] > '']
    #     con = sqlite3.connect('/home/mingyu/Bioinformatics/database/important_genes_from_Amlan.db')
    #     df.to_sql('ts', con, index=None, if_exists='replace')

    def to_sql(self):
        fpath = '/home/mingyu/Bioinformatics/database/important_genes_from_Amlan.txt'
        df = pd.read_csv(fpath, sep='\t', names=['miRNA', 'genes', 'corr'])
        genes = df['genes'].str.split(';')
        for idx in genes.index:
            print(idx + 1, genes.shape[0])
            if isinstance(genes[idx], float):
                continue

            gs = []
            for i, col in enumerate(genes[idx]):
                if ' /// ' in col:
                    gs.append(col.split(' /// ')[0])
                else:
                    gs.append(col)
            df.loc[idx, 'genes'] = ';'.join(gs)
        con = sqlite3.connect(fpath.replace('.txt', '.db'))
        df = df.sort_values('miRNA')
        df.to_sql('important_genes', con, index=None, if_exists='replace')


if __name__ == '__main__':
    hostname = socket.gethostname()
    if hostname == 'mingyu-Precision-Tower-7810':
        root = '/home/mingyu/Bioinformatics'
    elif hostname == 'DESKTOP-DLOOJR6' or hostname == 'DESKTOP-1NLOLK4':
        root = 'D:/Bioinformatics'
    elif hostname == 'mingyu-Inspiron-7559':
        root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
    else:
        root = '/lustre/fs0/home/mcha/Bioinformatics'

    rep = Report(root)
    # rep.get_compare_cell_tissue()
    # rep.stats_mannwhitneyu()
    rep.compare_others()
    # rep.missing_mir()