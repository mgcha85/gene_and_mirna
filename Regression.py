import pandas as pd
import numpy as np
import os
import socket
import pickle as pkl
from sklearn.linear_model import Lasso
from sklearn import linear_model

class Regression:
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
        self.scope = 500

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue_{}_score_vector_corr2.xlsx'.format(self.scope))
        df = pd.read_excel(fpath, index_col=0)

        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue_{}_score_vector.xlsx'.format(self.scope))
        df_gene = pd.read_excel(fpath, index_col=0)

        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'fantom_cage_by_tissue_mir_{}_score_vector.xlsx'.format(self.scope))
        df_mir = pd.read_excel(fpath, index_col=0)

        # contents = {}
        # for mir in df.index:
        #     genes = df.loc[mir, 'GENEs'].split(';')
        #     corr = df.loc[mir, 'corr (pearson)']
        #
        #     # mir_values = ';'.join(df_mir.loc[mir, :].round(4).astype(str))
        #     for gene in genes:
        #         if corr not in contents:
        #             contents[corr] = [df_gene.loc[gene, :].values]
        #         else:
        #             contents[corr].append(df_gene.loc[gene, :].values)
        #
        # with open('regression.cha', 'wb') as f:
        #     pkl.dump(contents, f)

        with open('regression.cha', 'rb') as f:
            contents = pkl.load(f)

        clf = linear_model.Lasso(alpha=0.1)
        for corr, genes in contents.items():
            corr = np.array(corr.split(';')).astype(float)
            clf.fit(genes, corr)
            Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
               normalize=False, positive=False, precompute=False, random_state=None,
               selection='cyclic', tol=0.0001, warm_start=False)
            print(clf.coef_)


if __name__ == '__main__':
    rg = Regression()
    rg.run()
