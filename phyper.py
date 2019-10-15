import sqlite3
import os
import pandas as pd
import socket
import numpy as np
from scipy.stats import hypergeom
from collections import OrderedDict


class phyper:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.n_gene = 4686
    #
    # def get_input_genes(self):
    #     fname = "data/human_cell_line_hCAGE_{}_score_vector.gct".format(500)
    #     df_gene = pd.read_csv(fname, sep="\t", skiprows=[0, 1], index_col=0)
    #     return df_gene.index
    #
    # def get_gene_set(self):
    #     dirname = '/home/mingyu/gsea_home/output/aug28/my_analysis.Gsea.1567011190321'
    #     flist = os.listdir(dirname)
    #     flist = [x for x in flist if x.endswith('.xls') and 'CHR' in x]
    #
    #     gene_set = []
    #     for i, fname in enumerate(flist):
    #         df_gsea = pd.read_csv(os.path.join(dirname, fname), sep='\t')
    #         gset = set(df_gsea['PROBE'].values)
    #         gene_set.append(set(gset))
    #     gene_set = set.union(*gene_set)
    #     return sorted(list(gene_set))

    def get_lasso_from_geneset(self):
        # GSEA
        dirname = '/home/mingyu/gsea_home/output/aug28/my_analysis.Gsea.1567011190321'
        flist = os.listdir(dirname)
        flist = [x for x in flist if x.endswith('.xls') and 'CHR' in x]

        gene_set = pd.Series()
        for fname in flist:
            fname_, ext = os.path.splitext(fname)
            df_gsea = pd.read_csv(os.path.join(dirname, fname), sep='\t')
            gene_set[fname_] = df_gsea['PROBE'].values

        # Lasso
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter.xlsx')
        df_las = pd.read_excel(fpath, index_col=0)

        gene_set_las = OrderedDict()
        for mir in df_las.index:
            genes = set(df_las.loc[mir, 'GENEs'].split(';'))
            coeffs = set(df_las.loc[mir, 'Coeffs'].split(';'))

            ser = pd.Series()
            for g, c in zip(genes, coeffs):
                ser[g] = c
            gene_set_las[mir] = ser

        # overlap between GSEA and Lasso
        df_res = pd.DataFrame(index=gene_set.index, columns=['GENEs', 'Comp', 'miRNA'])
        for sname in gene_set.index:
            gsea = gene_set[sname]
            gs_inter = []
            gs_comp = []
            for mir, gsl in gene_set_las.items():
                lasso = set(gsl.index)
                gsea = set(gsea)
                inter = set.intersection(lasso, gsea)

                if len(inter) > 1:
                    gs_inter.append(str(tuple(inter)))
                    gs_comp.append(str(tuple(lasso - inter)))

            if len(gs_inter) > 0:
                df_res.loc[sname, 'GENEs'] = ':'.join(gs_inter)
                df_res.loc[sname, 'Comp'] = ':'.join(gs_comp)
                df_res.loc[sname, 'miRNA'] = ':'.join(gene_set_las.keys())

        df_res = df_res.dropna(subset=['GENEs'])
        df_res.to_excel(fpath.replace('.xlsx', '_hyper.xlsx'))
    #
    # def get_lasso(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter.xlsx')
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     contents = []
    #     for i, mir in enumerate(df.index):
    #         genes = df.loc[mir, 'GENEs']
    #         genes = genes.split(';')
    #         contents.append(set(genes))
    #     return sorted(list(set.union(*contents)))
    #
    # def get_lasso_gsea(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter2.xlsx')
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     genes = pd.Series()
    #     for mir in df.index:
    #         gene_set = df.loc[mir, 'GENEs'].split(':')
    #         for gset in gene_set:
    #             gset = gset[1:-1].split(', ')
    #             gset = [x.strip().replace("'", '') for x in gset]
    #             genes[mir] = gset
    #     return genes
    #
    # def get_genes(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_500_score_corr2.xlsx')
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     contents = pd.Series()
    #     for i, mir in enumerate(df.index):
    #         genes = df.loc[mir, 'GENEs']
    #         genes = genes.split(';')
    #         contents[mir] = genes
    #     return contents
    #
    # def get_lasso_gene(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter.xlsx')
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     contents = pd.Series()
    #     for i, mir in enumerate(df.index):
    #         genes = df.loc[mir, 'GENEs']
    #         genes = genes.split(';')
    #         contents[mir] = genes
    #     return contents
    #
    # def intersection_gsea_lasso(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5', 'cell_lines/Regression_filter2.xlsx')
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     genes = []
    #     for mir in df.index:
    #         gene_set = df.loc[mir, 'GENEs'].split(':')
    #         for gset in gene_set:
    #             gset = gset[1:-1].split(', ')
    #             gset = [x.strip().replace("'", '') for x in gset]
    #             genes.append(set(gset))
    #     return sorted(list(set.union(*genes)))
    #
    # def get_gsea_lasso(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5', 'cell_lines/Regression_filter2.xlsx')
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     genes = pd.Series()
    #     for mir in df.index:
    #         gene_set = df.loc[mir, 'GENEs'].split(':')
    #         for gset in gene_set:
    #             gset = gset[1:-1].split(', ')
    #             gset = [x.strip().replace("'", '') for x in gset]
    #             genes[mir] = gset
    #     return genes
    #
    # def get_all_genes(self):
    #     fname = "data/human_cell_line_hCAGE_500_score_vector.gct".format()
    #     df_gene = pd.read_csv(fname, sep="\t", skiprows=[0, 1], index_col=0)
    #     return df_gene.index

    def get_param(self):
        dirname = '/home/mingyu/gsea_home/output/aug28/my_analysis.Gsea.1567011190321'
        flist = os.listdir(dirname)

        fmap = OrderedDict()
        for f in flist:
            if f.endswith('.xls') and 'CHR' in f:
                fmap[os.path.splitext(f)[0]] = os.path.join(dirname, f)

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter.xlsx')
        df_lasso = pd.read_excel(fpath, index_col=0)

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'Regression_filter_hyper.xlsx')
        df = pd.read_excel(fpath, index_col=0)

        df_res = pd.DataFrame(index=df.index)
        for fname in df.index:
            gset = df.loc[fname, 'GENEs'].split(':')
            comp = df.loc[fname, 'Comp'].split(':')
            mirna = df.loc[fname, 'miRNA'].split(':')

            df_gsea = pd.read_csv(os.path.join(dirname, fname + '.xls'), sep='\t')
            for genes, mir in zip(gset, mirna):
                genes = genes[1:-1].split(',')
                x = len(genes)
                m = df_gsea.shape[0]
                n = len(comp)
                k = len(df_lasso.loc[mir, 'GENEs'].split(';'))
                p = hypergeom.sf(x, n+m, m, k)
                df_res.loc[fname, mir] = ';'.join(map(str, [x, m, n, k, p]))

        df_res.dropna(how='all', axis=1).to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'p_values.xlsx'))

    def run(self):
        genes = self.get_genes()
        l_genes = self.get_lasso_gene()
        lg_genes = self.get_lasso_gsea()

        df = pd.concat([genes, l_genes, lg_genes], axis=1)
        df = df.dropna(how='any')
        df.columns = ['GENEs', 'L_GENEs', 'LG_GENEs']

        df_res = pd.DataFrame(data=np.zeros((df.shape[0], 5)), index=df.index, columns=['x', 'm', 'n', 'k', 'p-value'])
        for mir in df.index:
            x = len(df.loc[mir, 'LG_GENEs'])
            m = len(df.loc[mir, 'L_GENEs'])
            k = len(df.loc[mir, 'GENEs'])
            n = k - m
            df_res.loc[mir, 'x'] = x
            df_res.loc[mir, 'm'] = m
            df_res.loc[mir, 'k'] = k
            df_res.loc[mir, 'n'] = n
            df_res.loc[mir, 'p-value'] = hypergeom.sf(x, n+m, m, k)
        df_res.to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'p_values.xlsx'))


if __name__ == '__main__':
    p = phyper()
    p.get_param()
