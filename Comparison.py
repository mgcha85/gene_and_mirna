import sqlite3
import pandas as pd
import socket
import os
import numpy as np


class Comparison:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.bed_columns = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'rgb']
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    def processInput(self, df_ref, idx, fpath_rsc):
        if idx % 1000 == 0 or idx + 1 == df_ref.shape[0]:
            print('{:0.2f}%'.format(100 * (idx + 1) / df_ref.shape[0]))

        chromosome = df_ref.loc[idx, 'chromosome']
        strand = df_ref.loc[idx, 'strand']

        fpath, tname = fpath_rsc
        con = sqlite3.connect(fpath)
        if strand == '+':
            ref_tss = df_ref.loc[idx, 'start']
            rsc_tss = 'start'
        else:
            ref_tss = df_ref.loc[idx, 'end']
            rsc_tss = 'end'

        ref_start = ref_tss - 100
        ref_end = ref_tss + 100
        df_rsc = pd.read_sql_query("SELECT * FROM '{}_{}_{}' WHERE {} BETWEEN {} AND {}"
                                   "".format(tname, chromosome, strand, rsc_tss, ref_start, ref_end), con)

        if not df_rsc.empty:
            return idx, None
        else:
            return None, idx

    def run(self):
        from joblib import Parallel, delayed
        import multiprocessing

        ref_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.annotation.db')
        ref_con = sqlite3.connect(ref_path)
        df_ref = pd.read_sql_query("SELECT * FROM '{}' WHERE feature='transcript'".format('gencode.v30lift37'), ref_con)

        fpath_ens = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19.db')
        fpath_ucsc = os.path.join(self.root, 'database/UCSC/Genes', 'genes.db')
        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/hg19.cage_peak_phase1and2combined_coord.db')

        dfs_rsc = {'Ensembl': (fpath_ens, 'Ensembl'), 'UCSC': (fpath_ucsc, 'genes'), 'FANTOM': (fpath_fan, 'hg19_cage_peak_phase1and2combined_coord')}

        num_cores = multiprocessing.cpu_count() - 1
        for key, fpath_rsc in dfs_rsc.items():
            print(key)
            results = Parallel(n_jobs=num_cores)(delayed(self.processInput)(df_ref, idx, fpath_rsc) for idx in df_ref.index)

            corr_idx = []
            non_corr_idx = []
            for res in results:
                if res is None:
                    continue

                if res[0] is not None:
                    corr_idx.append(res[0])
                else:
                    non_corr_idx.append(res[1])
        #
        # for idx in df_ref.index:
        #     if idx % 1000 == 0 or idx + 1 == df_ref.shape[0]:
        #         print('{:0.2f}%'.format((100 * idx + 1) / df_ref.shape[0]))
        #
        #     chromosome = df_ref.loc[idx, 'chromosome']
        #     strand = df_ref.loc[idx, 'strand']
        #
        #     if strand == '+':
        #         ref_tss = df_ref.loc[idx, 'start']
        #     else:
        #         ref_tss = df_ref.loc[idx, 'end']
        #
        #     if chromosome not in df_ens_grp:
        #         continue
        #     df_ens_chr = df_ens_grp.get_group(chromosome)
        #     df_ens_chr_str = df_ens_chr[df_ens_chr['Strand'] == strand]
        #
        #     if strand == '+':
        #         src_tss = df_ens_chr_str['Transcript start (bp)']
        #     else:
        #         src_tss = df_ens_chr_str['Transcript end (bp)']
        #
        #     df_ens_chr_str_corr = df_ens_chr_str[(src_tss - ref_tss).abs() < 100]
        #     if not df_ens_chr_str_corr.empty:
        #         corr_idx.append(idx)
        #     else:
        #         non_corr_idx.append(idx)

            df_ref.loc[corr_idx, '<100bp ({})'.format(key)] = 'O'
            df_ref.loc[non_corr_idx, '<100bp ({})'.format(key)] = 'X'

        out_path = os.path.join('{}.xlsx'.format(self.__class__.__name__))
        writer = pd.ExcelWriter(out_path, engine='xlsxwriter')
        df_ref.to_excel(writer, sheet_name='corresponding', index=None)
        writer.save()
        writer.close()

    def plot_distance_histogram(self):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database', 'tss_comparison.xlsx')
        df = pd.read_excel(fpath)
        df_fan = df[df['distance_fan'] < 1000]
        df_ucsc = df[df['distance_ucsc'] < 1000]

        print('FANTOM & Ensembl: {:0.2f}% [<1kb]'.format(100 * df_fan.shape[0] / df.shape[0]))
        print('UCSC & Ensembl: {:0.2f}% [<1kb]'.format(100 * df_ucsc.shape[0] / df.shape[0]))

        bins = [1]
        for i in range(2, 6):
            bins.append(10 ** i)
        bins.append(10 ** 10)

        hist_fan, bin_fan = np.histogram(df['distance_fan'].values, bins=bins)
        hist_ucsc, bin_ucsc = np.histogram(df['distance_ucsc'].values, bins=bins)

        plt.subplot(211)
        xaxis = np.log10(bins[:-1])
        plt.bar(xaxis, hist_fan / df.shape[0])
        plt.xticks(xaxis, bins[:-1])
        plt.title('Ensembl vs. Fantom')
        plt.grid()
        plt.legend()

        plt.subplot(212)
        plt.bar(xaxis, hist_ucsc / df.shape[0])
        plt.xticks(xaxis, bins[:-1])
        plt.title('Ensembl vs. UCSC')
        plt.grid()
        plt.show()

    def split_table(self):
        # fpath_ens = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19.db')
        fpath_ucsc = os.path.join(self.root, 'database/UCSC/Genes', 'genes.gtf')
        # fpath_fan = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.db')
        df = pd.read_csv(fpath_ucsc, sep='\t', names=self.gtf_columns, comment='#')
        con = sqlite3.connect(fpath_ucsc.replace('.gtf', '.db'))
        df.to_sql('UCSC', con)
        exit(1)
        con = sqlite3.connect(fpath_fan)
        df = pd.read_sql_query("SELECT * FROM 'Fantom'", con)

        tname = 'Fantom'
        con_out = sqlite3.connect(fpath_fan.replace('.db', '_out.db'))
        # con = sqlite3.connect(fpath_ucsc)
        # df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
        df_chr = df.groupby('chromosome')
        for chr, df_sub in df_chr:
            if len(chr) > 5:
                continue
            for strand, df_sub_sub in df_sub.groupby('strand'):
                df_sub_sub.to_sql('{}_{}_{}'.format(tname, chr, strand), con_out, index=None)


if __name__ == '__main__':
    comp = Comparison()
    comp.split_table()
