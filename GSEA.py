import pandas as pd
import sqlite3
import os
import socket
import gseapy as gp
import numpy as np


class GSEA:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.scope = 500

    def filter_by_hgnc(self):
        fname = "data/human_cell_line_hCAGE_{}_score_vector.gct".format(self.scope)
        df_gene = pd.read_csv(fname, sep="\t", skiprows=[0, 1], index_col=0)
        df_hgnc = pd.read_csv("data/hgnc-symbol-check.csv", skiprows=[0], index_col=0)
        df_hgnc = df_hgnc.dropna(subset=['Approved symbol'])
        df_hgnc = df_hgnc.loc[~df_hgnc.index.duplicated(keep='first')]

        for gene in df_gene.index:
            if gene in df_hgnc.index:
                symbol = df_hgnc.loc[gene, 'Approved symbol']
                df_gene.loc[gene, 'NAME'] = symbol
            else:
                df_gene = df_gene.drop(gene)

        df_gene = df_gene.set_index('NAME', drop=True)
        df_gene = df_gene.reset_index()
        contents = ["#1.2", "{} {}".format(df_gene.shape[0], df_gene.shape[1] - 2), '\t'.join(df_gene.columns)]
        for idx in df_gene.index:
            contents.append('\t'.join(df_gene.loc[idx].astype(str)))
        with open(fname, 'wt') as f:
            f.write('\n'.join(contents))

    def set_data(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector.xlsx'.format(self.scope))
        df = pd.read_excel(fpath)
        df = df.dropna(how='any', axis=0)
        df.loc[:, 'DESCRIPTION'] = np.nan
        df = df.rename(columns={'gene': 'NAME'})

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_mir_{}_score_vector.xlsx'.format(self.scope))
        df_mir = pd.read_excel(fpath)
        df_mir = df_mir.dropna(how='any', axis=0)
        df_mir.loc[:, 'DESCRIPTION'] = np.nan
        df_mir = df_mir.rename(columns={'gene': 'NAME'})

        df = pd.concat([df, df_mir]).reset_index(drop=True)

        columns = ['NAME', 'DESCRIPTION'] + list(df.columns)[1:-1]
        df = df[columns]
        df.iloc[:, 2:] = df.iloc[:, 2:].round(2)

        contents = ["#1.2", "{} {}".format(df.shape[0], df.shape[1] - 2), '\t'.join(df.columns)]
        for idx in df.index:
            contents.append('\t'.join(df.loc[idx].astype(str)))
        with open('data/human_cell_line_hCAGE_{}_score_vector.gct'.format(self.scope), 'wt') as f:
            f.write('\n'.join(contents))

    def get_class(self):
        df = pd.read_csv('data/human_cell_line_hCAGE_{}_score_vector.gct'.format(self.scope), sep='\t', skiprows=[0, 1])
        cls = df.iloc[:, 2:].mean(axis=0).astype(int).astype(str)
        # cls[-N // 2:] = 1

        contents = ['#numeric', '#PeakProfle', ' '.join(cls)]
        with open('data/human_cell_line_hCAGE_{}_score_vector.cls'.format(self.scope), 'wt') as f:
            f.write('\n'.join(contents))

    # def get_class(self):
    #     df = pd.read_csv('data/human_cell_line_hCAGE_100_score_vector.gct', sep='\t', skiprows=[0, 1])
    #     N = df.shape[1] - 2
    #
    #     cls = np.zeros(N).astype(int)
    #     cls[-N // 2:] = 1
    #     cls = cls.astype(str)
    #
    #     contents = ['{} 2 1'.format(N), '#MUT WT', ' '.join(cls)]
    #     with open('data/human_cell_line_hCAGE_100_score_vector.cls', 'wt') as f:
    #         f.write('\n'.join(contents))

    def test(self):
        df_p53 = pd.read_csv("data/P53_collapsed_symbols.gct", sep="\t", skiprows=[0, 1])
        df_my = pd.read_csv("data/human_cell_line_hCAGE_100_score_vector.gct", sep="\t", skiprows=[0, 1])
        df_my['NAME'] = df_p53['NAME'][:df_my.shape[0]]
        df_my.to_csv("data/human_cell_line_hCAGE_100_score_vector.gct", sep="\t", index=None)

    def get_mir_gene(self):
        dirname = '/home/mingyu/gsea_home/output/jul31/my_analysis.Gsea.1564597163164'
        flist = [x for x in os.listdir(dirname) if x.endswith('.xls') and 'CHR' in x]
        for fname in flist:
            fpath = os.path.join(dirname, fname)
            df = pd.read_csv(fpath, sep='\t')
            df_mir = df[df['PROBE'].str.contains('MIR')]
            if not df_mir.empty:
                print(fname)

    def compare_amlan(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter2.xlsx')
        df_res = pd.read_excel(fpath, index_col=0)

        fpath_ref = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'important_genes2.xlsx')
        df_reg = pd.read_excel(fpath_ref, index_col=0)
        df_reg_grp = df_reg.groupby('miRNA')

        for mir in df_res.index:
            gene_set = df_res.loc[mir, 'GENEs'].split(':')
            if mir in df_reg_grp.groups:
                df_reg_pos = df_reg_grp.get_group(mir)
                df_reg_pos = df_reg_pos.set_index('gene')

                gene_inter = []
                weights = []
                for gs in gene_set:
                    gs = gs[1:-1].split(', ')
                    gs = [x.strip().replace("'", '') for x in gs]

                    inter = list(set.intersection(set(gs), set(df_reg_pos.index)))
                    gene_inter.append(str(tuple(inter)))
                    weights.append(str(tuple(df_reg_pos.loc[inter, 'weight'].values.astype(str))))

                gene_inter = ':'.join(gene_inter)
                weights = ':'.join(weights)
                df_res.loc[mir, 'Truth'] = gene_inter
                df_res.loc[mir, 'Coeff (Truth)'] = weights
        df_res.to_excel(fpath)

    def lasso_gsea(self):
        # Lasso
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter.xlsx')
        df_las = pd.read_excel(fpath, index_col=0)

        gene_set_las = {}
        for mir in df_las.index:
            genes = set(df_las.loc[mir, 'GENEs'].split(';'))
            coeffs = set(df_las.loc[mir, 'Coeffs'].split(';'))

            ser = pd.Series()
            for g, c in zip(genes, coeffs):
                ser[g] = c
            gene_set_las[mir] = ser

        # GSEA
        dirname = '/home/mingyu/gsea_home/output/aug28/my_analysis.Gsea.1567011190321'
        flist = os.listdir(dirname)
        flist = [x for x in flist if x.endswith('.xls') and 'CHR' in x]

        gene_set = []
        for fname in flist:
            df_gsea = pd.read_csv(os.path.join(dirname, fname), sep='\t')
            gene_set.append(set(df_gsea['PROBE'].values))

        # overlap between GSEA and Lasso
        df_res = pd.DataFrame(index=gene_set_las.keys(), columns=['GENEs', 'Coeffs'])
        for mir, gsl in gene_set_las.items():
            # gs_inter = set()
            gs_inter = []
            gs_coeff = []
            for gs in gene_set:
                inter = set.intersection(set(gsl.index), gs)
                if len(inter) > 1:
                    gs_inter.append(str(tuple(inter)))
                    gs_coeff.append(str(tuple(gsl.loc[inter])))
                # gs_inter = set.union(gs_inter, set.intersection(gsl, gs))

            if len(gs_inter) > 0:
                df_res.loc[mir, 'GENEs'] = ':'.join(gs_inter)
                df_res.loc[mir, 'Coeffs'] = ':'.join(gs_coeff)
        df_res = df_res.dropna(subset=['GENEs'])
        df_res.to_excel(fpath.replace('.xlsx', '2.xlsx'))

    def chart_plot(self):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter2.xlsx')
        df = pd.read_excel(fpath, index_col=0)

        buffer = np.zeros(df.shape[0])
        xaxis = np.arange(df.shape[0])

        ser = pd.Series()
        for i, mir in enumerate(df.index):
            genes = df.loc[mir, 'GENEs'].split(';')
            coeffs_ = np.array(list(map(float, df.loc[mir, 'Coeffs'].split(';'))))
            buffer[i] = len(genes)
            mean = coeffs_.mean()
            std = coeffs_.std()

            for c, g in zip(coeffs_, genes):
                ser[';'.join([mir, g])] = c

        plt.subplot(211)
        plt.bar(xaxis, buffer)
        plt.title('[COEFFS] mean: {:0.2f}, std: {:0.2f}'.format(mean, std))
        plt.plot([xaxis[0], xaxis[-1]], mean * np.ones(2), color='r', linestyle=':')
        plt.xticks(xaxis[::10], df.index[::10], rotation=30, fontsize=6)
        plt.grid()

        df = df.dropna(subset=['Truth'])
        buffer = np.zeros(df.shape[0])
        xaxis = np.arange(df.shape[0])
        for i, mir in enumerate(df.index):
            genes = df.loc[mir, 'Truth'].split(';')
            coeffs = np.array(list(map(float, df.loc[mir, 'Coeff (Truth)'].split(';'))))
            buffer[i] = len(genes)
            # coeffs = []
            # for g in genes:
            #     coeffs.append(ser[';'.join([mir, g])])
            coeffs = np.array(coeffs)
            mean = coeffs.mean()
            std = coeffs.std()

        plt.subplot(212)
        plt.bar(xaxis, buffer)
        plt.title('[COEFFS] mean: {:0.2f}, std: {:0.2f}'.format(mean, std))
        plt.plot([xaxis[0], xaxis[-1]], mean * np.ones(2), color='r', linestyle=':')
        plt.xticks(xaxis, df.index, rotation=30, fontsize=6)
        plt.grid()
        plt.show()

    def run(self):
        phenoA, phenoB, class_vector = gp.parser.gsea_cls_parser("data/P53.cls")
        gene_exp = pd.read_csv("data/human_cell_line_hCAGE_100_score_vector.gct", sep="\t", skiprows=[0, 1])

        # gene_exp = gene_exp[['NAME', 'DESCRIPTION', *gene_exp.columns[1:-1]]]
        # gene_exp.to_csv("data/human_cell_line_hCAGE_100_score_vector.txt", sep="\t", index=None)

        gs_res = gp.gsea(data=gene_exp,  # or data='./P53_resampling_data.txt'
                         gene_sets='KEGG_2016',  # enrichr library names
                         cls='data/human_cell_line_hCAGE_100_score_vector.cls',  # cls=class_vector
                         # set permutation_type to phenotype if samples >=15
                         permutation_type='phenotype',
                         permutation_num=100,  # reduce number to speed up test
                         outdir='data',
                         no_plot=True,  # Skip plotting
                         method='signal_to_noise',
                         processes=4,
                         format='png')
        print(gs_res.res2d.head())


if __name__ == '__main__':
    # import urllib.request
    # url = 'https://raw.githubusercontent.com/stephens999/ash-orig/master/data/GSEA/p53/P53_collapsed_symbols.gct'
    # urllib.request.urlretrieve(url, 'P53_collapsed_symbols.gct')
    gs = GSEA()
    # gs.set_data()
    # gs.filter_by_hgnc()
    # gs.get_class()
    gs.lasso_gsea()
    # gs.compare_amlan()
    # gs.chart_plot()
