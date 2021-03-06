import os
import pandas as pd
import numpy as np
import sqlite3
import socket
from Database import Database
import sys
from Server import Server
import zipfile


class data_preparation:
    def __init__(self, root):
        self.root = root

    def compress(self, in_file, out_file):
        compression = zipfile.ZIP_DEFLATED
        zf = zipfile.ZipFile(out_file, mode="w")
        try:
            print("compress '{}'...".format(in_file))
            zf.write(in_file, os.path.split(in_file[1]), compress_type=compression)
        except FileNotFoundError:
            print("An error occurred")
        finally:
            zf.close()

    def to_server(self, root, tag):
        from Server import Server
        import sys

        which = 'newton'
        server = Server(root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname__ = os.path.split(local_path)
        fname = fname__.replace('.py', '_{}.py'.format(tag))

        curdir = os.getcwd().split(os.sep)[-1]
        server_root = os.path.join(server.server, 'source', curdir).replace(os.sep, '/')
        server_path = local_path.replace(dirname, server_root)
        server_path = server_path.replace(fname__, fname)

        server.job_script(fname, src_root=server_root, time='08:00:00', pversion=3)
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', server_root + '/' + 'dl-submit.slurm')

        stdin, stdout, stderr = server.ssh.exec_command("cd {root};sbatch {root}/dl-submit.slurm".format(root=server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    # 1. data download from FANTOM5
    def download_cell_line_fantom(self):
        import urllib.request

        df = pd.read_csv('fantom_cell_line_list.txt')
        for idx in df.index:
            link = df.loc[idx, 'link']
            url, fname = os.path.split(link)
            dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
            fpath = os.path.join(dirname, fname)
            urllib.request.urlretrieve(link, fpath)
            # df.loc[idx, 'link'] = link.replace('%', '%25')
        # df.to_csv('fantom_cell_line_list.txt', index=None)

    # 2. check if it is hg19 or not
    def check_chain(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        fpath = os.path.join(dirname, 'list/00_human.cell_line.hCAGE.hg19.assay_sdrf2.xlsx')
        df_list = pd.read_excel(fpath)
        df_sub = df_list[~df_list['File Name.1'].str.contains('hg19')]
        print(df_sub)

    def merge_processInput(self, row, columns):
        fpath, cline, ename = row
        return pd.read_csv(fpath, sep='\t', compression='gzip', names=columns), cline, ename

    # 3. integrate cell line data
    def merge_cline_db(self):
        from joblib import Parallel, delayed
        num_cores = 6

        def as_batch(data, i, batch_size=100):
            N = len(data)
            sidx = i * batch_size
            eidx = sidx + batch_size
            if eidx > N:
                eidx = N
            return data[sidx: eidx]

        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        fpath = os.path.join(dirname, 'list/00_human.cell_line.hCAGE.hg19.assay_sdrf2.xlsx')
        df_list = pd.read_excel(fpath)
        df_grp = df_list.groupby('cell line')
        columns = ['chromosome', 'start', 'end', 'location', 'score', 'strand']

        con = sqlite3.connect(os.path.join(dirname, 'human_cell_line_hCAGE.db'))
        fpaths = []
        fails = []
        multi = []
        for cline, df_sub in df_grp:
            if df_sub.shape[0] > 1:
                print(cline, df_sub.shape[0])
                multi.append([cline, df_sub.shape[0]])
            idx = df_sub.index[0]
            fname = df_sub.loc[idx, 'File Name.1']
            ename = df_sub.loc[idx, 'Extract Name']
            fname = fname.replace('.nobarcode.bam', '.ctss.bed.gz')

            fpath = os.path.join(dirname, fname)
            if not os.path.exists(fpath):
                fails.append(fname)
                continue
            fpaths.append([fpath, cline, ename])

        with open('fails.txt', 'wt') as f:
            f.write('\n'.join(fails))
        df_multi = pd.DataFrame(multi, columns=['cell line', 'N'])
        df_multi.to_csv('multiple_cell_lines.csv', index=None)

        N = int((len(fpaths) + num_cores) // num_cores)
        for i in range(N):
            dfs = Parallel(n_jobs=num_cores)(delayed(self.merge_processInput)(row, columns) for row in as_batch(fpaths, i, num_cores))
            for row in dfs:
                df, cline, ename = row
                df[['chromosome', 'start', 'end', 'strand', 'score']].to_sql(cline, con, if_exists='replace', index=None)

    def cage_tags_to_db(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        fpath = os.path.join(dirname, 'robust_phase1_pls_2.tpm.desc121113.osc.txt.gz.tmp')

        # get comment
        with open(fpath, 'rt') as f:
            lines = f.read(1<<20).split('\n')

        comments = []
        for line in lines:
            if line[:2] == '##':
                comments.append(line)
            else:
                print('enough')
                break

        with open(os.path.join(dirname, 'comments.txt'), 'wt') as f:
            f.write('\n'.join(comments))

        # get contents
        con = sqlite3.connect(os.path.join(dirname, 'robust_phase1_pls_2.tpm.desc121113.osc.txt.gz.db'))
        for i, df_sub in enumerate(pd.read_csv(fpath, sep='\t', iterator=True, chunksize=1<<16, comment='#')):
            print('batch {}'.format(i))
            df_sub.to_sql('overall', con, if_exists='append')

    def split_cage_tags(self):
        root = os.path.join(self.root, 'database/Fantom/v5/tissues')
        dirList = [x for x in os.listdir(root) if os.path.isdir(os.path.join(root, x)) and x != 'out']

        result = []
        for dir in dirList:
            dirname = os.path.join(root, dir)
            flist = [x for x in os.listdir(dirname) if x.endswith('.gz')]
            for fname in flist:
                result.append([dir, fname.split('.')[2]])
        df_list = pd.DataFrame(result, columns=['tissue', 'CAGE FF_ID'])

        fpath = os.path.join(root, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz')
        out_path = os.path.join(root, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.db')
        con = sqlite3.connect(out_path)
        for df_sub in pd.read_csv(fpath, sep='\t', iterator=True, chunksize=1 << 16, comment='#'):
            df_sub = df_sub.dropna(subset=['description'])
            df_grp = df_list.groupby('tissue')
            print(len(df_grp))

            for tissue__, df_tis in df_grp:
                columns__ = []
                for idx in df_tis.index:
                    cid = df_tis.loc[idx, "CAGE FF_ID"]
                    columns__ += [col for col in df_sub.columns[7:] if cid in col]
                columns__ = list(set(columns__))
                if columns__:
                    columns = list(df_sub.columns[:7]) + columns__
                    df_sub[columns].to_sql(tissue__, con, if_exists='append', index=None)
                else:
                    print('{} is not available'.format(tissue__))
        self.compress(out_path, out_path.replace('.db', '.zip'))

    # def split_cage_tags(self):
    #     dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
    #
    #     fpath = os.path.join(dirname, 'tissues_fantom_rna.xlsx')
    #     df_list = pd.read_excel(fpath)
    #
    #     df_ref = pd.read_excel(os.path.join(self.root, "Papers/Lasso", "supp_gkv608_nar-00656-h-2015-File009_2.xls"))
    #     df_ref = df_ref.set_index("tissue")
    #
    #     fpath = os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz')
    #     con = sqlite3.connect(os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.db'))
    #     for df_sub in pd.read_csv(fpath, sep='\t', iterator=True, chunksize=1 << 10, comment='#'):
    #         df_sub = df_sub.dropna(subset=['description'])
    #         for tissue__ in df_list['FANTOM']:
    #             # tissue = tissue__.replace('_', '%20')
    #             cage_id = df_ref.loc[tissue__, "CAGE FF_ID"]
    #             if isinstance(cage_id, pd.Series):
    #                 columns__ = []
    #                 for cid in cage_id:
    #                     columns__ += [col for col in df_sub.columns[7:] if cid in col]
    #                 columns__ = list(set(columns__))
    #             else:
    #                 columns__ = [col for col in df_sub.columns[7:] if cage_id in col]
    #             if columns__:
    #                 columns = list(df_sub.columns[:7]) + columns__
    #                 df_sub[columns].to_sql(tissue__, con, if_exists='append', index=None)
    #             else:
    #                 print('{} is not available'.format(tissue__))

    def split_cage_tags2(self, download=False):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')

        if download:
            import urllib.request
            url = 'https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz'

            _, fname = os.path.split(url)
            urllib.request.urlretrieve(url, os.path.join(dirname, fname))

        fpath = os.path.join(dirname, 'list', '00_human.cell_line.hCAGE.hg19.assay_sdrf2.xlsx')
        df_list = pd.read_excel(fpath)
        cell_lines = set(df_list['cell line'].dropna())

        no_data = set()
        fpath = os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz')
        con = sqlite3.connect(os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.db'))
        for df_sub in pd.read_csv(fpath, sep='\t', iterator=True, chunksize=1 << 16, comment='#'):
            df_sub = df_sub.dropna(subset=['description'])
            for cline in cell_lines:
                if cline in no_data:
                    continue

                cline = cline.replace(' ', '%20')
                cline = cline.replace('/', '%2f')
                cline = cline.replace('.', '%2e')
                cline = cline.replace('^', '%5e')
                cline = cline.replace('(', '%28')
                cline = cline.replace(')', '%29')
                columns__ = [col for col in df_sub.columns[7:] if cline in col]
                if columns__:
                    columns__ = sorted(list(set(columns__)))
                    columns = list(df_sub.columns[:7]) + columns__
                    df_sub[columns].to_sql(cline, con, if_exists='append', index=None)
                else:
                    print('{} is not available'.format(cline))
                    if cline not in no_data:
                        no_data.add(cline)

        with open('not_avail_cline', 'wt') as f:
            f.write('\n'.join(sorted(list(no_data))))

    def avg_cage_tags(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        con = sqlite3.connect(os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.db'))

        out_path = os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc_avg.db')
        con_out = sqlite3.connect(out_path)
        for tissue in Database.load_tableList(con):
            df = pd.read_sql("SELECT * FROM '{}'".format(tissue), con)
            df['avg_score'] = df.iloc[:, 7:].mean(axis=1)

            location1 = df[df.columns[0]].str.split(":", expand=True)
            location2 = location1[1].str.split(",", expand=True)
            location3 = location2[0].str.split(".", expand=True)
            df['chromosome'] = location1[0]
            df['strand'] = location2[1]
            df['start'] = location3[0].astype(int)
            df['end'] = location3[2].astype(int)

            columns = ['chromosome', 'start', 'end', 'strand'] + list(df.columns[1:7]) + ['avg_score']
            df[columns].to_sql(tissue, con_out, if_exists='replace', index=None)
        self.compress(out_path, out_path.replace('.db', '.zip'))

    def avg_cage_tags2(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        con = sqlite3.connect(os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.db'))
        con_out = sqlite3.connect(os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc_avg.db'))
        for tissue in Database.load_tableList(con):
            df = pd.read_sql("SELECT * FROM '{}'".format(tissue), con)
            df['avg_score'] = df.iloc[:, 7:].mean(axis=1)

            location1 = df[df.columns[0]].str.split(":", expand=True)
            location2 = location1[1].str.split(",", expand=True)
            location3 = location2[0].str.split(".", expand=True)
            df['chromosome'] = location1[0]
            df['strand'] = location2[1]
            df['start'] = location3[0].astype(int)
            df['end'] = location3[2].astype(int)

            columns = ['chromosome', 'start', 'end', 'strand'] + list(df.columns[1:7]) + ['avg_score']
            df[columns].to_sql(tissue, con_out, if_exists='replace', index=None)

    def set_cage_data(self):
        # tissue list
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'tissues_fantom_rna.xlsx')
        df_tis = pd.read_excel(fpath, sheet_name='Sheet1', index_col=1)

        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        con = sqlite3.connect(os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc_avg.db'))
        con_out = sqlite3.connect(os.path.join(dirname, 'CAGE_tag_tissue_spt.db'))
        for tissue in Database.load_tableList(con):
            tissue_rna = df_tis.loc[tissue, 'RNA-seq']
            df = pd.read_sql("SELECT chromosome, start, end, strand, avg_score as score FROM '{}' ORDER BY chromosome, start".format(tissue), con)
            for str, df_str in df.groupby('strand'):
                for chr, df_chr in df_str.groupby('chromosome'):
                    df_chr.drop(['chromosome', 'strand'], axis=1).to_sql('{}_{}_{}'.format(tissue_rna, chr, str), con_out, if_exists='replace', index=None)

    def set_gencode(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v32lift37.annotation_attr.db')
        con = sqlite3.connect(fpath)
        con_out = sqlite3.connect(fpath.replace('.db', '_spt.db'))
        for tname in Database.load_tableList(con):
            chr, str = tname.split('_')
            df = pd.read_sql("SELECT start, end, gene_name, transcript_id FROM '{}' WHERE feature='transcript' AND "
                             "gene_type='protein_coding' AND transcript_type='protein_coding'".format(tname), con)
            if str == '+':
                df['tss'] = df['start']
                df['start'] = df['tss'] - 100
                df['end'] = df['tss'] + 100
            else:
                df['tss'] = df['end']
                df['start'] = df['tss'] - 100
                df['end'] = df['tss'] + 100
            df.drop('tss', axis=1).to_sql(tname, con_out, if_exists='replace', index=None)

    def download_tissue_fantom(self):
        import urllib.request

        url__ = 'http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.tissue.hCAGE/{}'
        dirname = os.path.join(self.root, 'database/Fantom/v5/tissues')
        fpath = os.path.join(dirname, '00_human.tissue.hCAGE.hg19.assay_sdrf.csv')
        df = pd.read_csv(fpath, sep='\t')
        df['Comment [sample_name]'] = df['Comment [sample_name]'].str.split(', ').str[0]
        df['Comment [sample_name]'] = df['Comment [sample_name]'].str.split(' - ').str[0]

        for grp, df_grp in df.groupby('Comment [sample_name]'):
            subdir = os.path.join(dirname, grp)
            if not os.path.exists(subdir):
                os.mkdir(subdir)
            for idx in df_grp.index:
                fname = df_grp.loc[idx, 'File Name.1'].replace('.nobarcode.bam', '.ctss.bed.gz').replace('%', '%25')
                print(grp, fname)
                url = url__.format(fname)
                urllib.request.urlretrieve(url, os.path.join(subdir, fname))

    def correlation_replicates(self, df, cols):
        from itertools import combinations

        report = pd.Series()
        for comb in combinations(cols, 2):
            coeff = np.corrcoef(df[comb[0]].fillna(0), df[comb[1]].fillna(0))
            report['&'.join(comb)] = coeff[0, 1]
        return report

    def avg_fantom_by_tissue(self):
        columns = ['chromosome', 'start', 'end', 'loc', 'score', 'strand']
        dirname_ = os.path.join(self.root, 'database/Fantom/v5/tissues')

        fpath = os.path.join(dirname_, 'FANTOM_tissue.db')
        con_out = sqlite3.connect(fpath)

        # writer = pd.ExcelWriter(os.path.join(self.root, 'database/Fantom/v5/tissues', 'duplicates.xlsx'), engine='xlsxwriter')
        flist = [os.path.join(dirname_, x) for x in os.listdir(dirname_)]
        dirlist = sorted([x for x in flist if os.path.isdir(x)])

        for i, dirname in enumerate(dirlist):
            tissue = dirname.split('/')[-1]
            dfs = []
            index = []

            flist = os.listdir(dirname)
            flist = [f for f in flist if f.endswith('.gz')]
            print(i, dirname)

            for fname in flist:
                fpath = os.path.join(dirname, fname)
                df_fid = pd.read_csv(fpath, compression='gzip', sep='\t', names=columns)
                df_fid.index = df_fid['chromosome'] + ':' + df_fid['start'].astype(str) + '-' + df_fid['end'].astype(str)
                index.append(set(df_fid.index))
                dfs.append(df_fid)

            index = sorted(list(set.union(*index)))
            df_res = pd.DataFrame(index=index, columns=df_fid.columns)
            cols = []
            for i, df_fid in enumerate(dfs):
                cols.append('score (Rep {})'.format(i + 1))
                df_res.loc[df_fid.index, cols[-1]] = df_fid['score'].astype(float)
                df_res.loc[df_fid.index, cols[-1]] /= df_res.loc[df_fid.index, cols[-1]].sum() / 1e6
                df_res.loc[df_fid.index, columns] = df_fid[columns]

            # report = self.correlation_replicates(df_res, cols)
            # report.to_excel(writer, sheet_name=tissue)
            df_res['score'] = df_res[cols].mean(axis=1)
            print(df_res[cols].sum(axis=0))
            df_res.drop(['loc'], axis=1).to_sql(tissue, con_out, if_exists='replace', index=None)
        # writer.save()
        # writer.close()

    def rna_cage(self):
        fpath = os.path.join(self.root, 'Papers/Lasso', 'supp_gkv608_nar-00656-h-2015-File009.xls')
        df = pd.read_excel(fpath, skiprows=1)
        for idx in df.index:
            tissue = df.loc[idx, 'RNASEQ samples'].split('.')[1].split('_')[0]
            df.loc[idx, 'tissue'] = tissue

        dfs = []
        for tis, df_sub in df.groupby('tissue'):
            df_sub[["CAGE FF_ID", "CAGE samples"]] = df_sub[["CAGE FF_ID", "CAGE samples"]].fillna(method='ffill')
            dfs.append(df_sub)
        pd.concat(dfs).to_excel(fpath.replace('.xls', '_2.xls'), index=None)


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

    dp = data_preparation(root)
    if hostname == 'DESKTOP-DLOOJR6' or hostname == 'DESKTOP-1NLOLK4':
        dp.to_server(root, "")
    else:
        dp.split_cage_tags()
        # dp.avg_cage_tags()
        # dp.set_cage_data()

        # dp.split_cage_tags2()
        # dp.avg_cage_tags2()
        # dp.rna_cage()
        # exit(1)
        # dp.split_cage_tags()
        # dp.avg_cage_tags()
        # dp.set_gencode()

        # dp.avg_fantom_by_tissue()

        # from Correlation_gpu import Correlation
        # from Regression import Regression
        #
        # cor = Correlation()
        # rg = Regression()
        #
        # cor.sum_fan(100, ref='gene')
        # cor.sum_fan(100, ref='mir')
        # rg.regression(100, 'nz')
