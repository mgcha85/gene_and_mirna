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

    def run(self):
        fpath = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19.db')
        con = sqlite3.connect(fpath)
        df_ens = pd.read_sql_query("SELECT * FROM 'Ensembl'", con)

        fpath = os.path.join(self.root, 'database/UCSC/Genes', 'genes.db')
        con = sqlite3.connect(fpath)
        df_ucsc = pd.read_sql_query("SELECT * FROM protein_coding", con)
        df_ucsc_grp = df_ucsc.groupby('seqname')

        fpath = os.path.join(self.root, 'database/Fantom/v5/hg19.cage_peak_phase1and2combined_coord.bed')
        df_fan = pd.read_csv(fpath, sep='\t', names=self.bed_columns)
        df_fan_grp = df_fan.groupby('chromosome')

        contents = []
        for idx in df_ens.index:
            if idx % 1000 == 0 or idx + 1 == df_ens.shape[0]:
                print('{:0.2f}%'.format((100 * idx + 1) / df_ens.shape[0]))

            start = df_ens.loc[idx, 'Transcript start (bp)']
            end = df_ens.loc[idx, 'Transcript end (bp)']
            if df_ens.loc[idx, 'Strand'] < 0:
                strand = '-'
            else:
                strand = '+'
            tname = df_ens.loc[idx, 'Transcript name']

            chrom = df_ens.loc[idx, 'Chromosome/scaffold name']
            if len(chrom) > 5:
                continue

            df_fan_chr = df_fan_grp.get_group(chrom)
            df_fan_chr_str = df_fan_chr[df_fan_chr['strand'] == strand]

            df_ucsc_chr = df_ucsc_grp.get_group(chrom)
            df_ucsc_chr_str = df_ucsc_chr[df_ucsc_chr['strand'] == strand]

            if strand == '+':
                distance_fan = (df_fan_chr_str['start'] - start).abs()
                midx_fan = distance_fan.idxmin()
                closest_fan = df_fan_chr_str.loc[midx_fan, 'start']

                distance_ucsc = (df_ucsc_chr_str['start'] - start).abs()
                midx_ucsc = distance_ucsc.idxmin()
                closest_ucsc = df_ucsc_chr_str.loc[midx_ucsc, 'start']

            else:
                distance_fan = (df_fan_chr_str['end'] - end).abs()
                midx_fan = distance_fan.idxmin()
                closest_fan = df_fan_chr_str.loc[midx_fan, 'end']

                distance_ucsc = (df_ucsc_chr_str['end'] - end).abs()
                midx_ucsc = distance_ucsc.idxmin()
                closest_ucsc = df_ucsc_chr_str.loc[midx_ucsc, 'end']

            contents.append([chrom, start, end, strand, tname, closest_fan, distance_fan[midx_fan], closest_ucsc, distance_ucsc[midx_ucsc]])
        df_res = pd.DataFrame(contents, columns=['chrom', 'start', 'end', 'strand', 'tname', 'closest_fan', 'distance_fan', 'closest_ucsc', 'distance_ucsc'])
        df_res.to_excel(os.path.join(self.root, 'database', 'tss_comparison.xlsx'), index=None)

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


if __name__ == '__main__':
    comp = Comparison()
    comp.plot_distance_histogram()
