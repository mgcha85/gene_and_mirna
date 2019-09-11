import pandas as pd
import numpy as np
import os
import socket
import pickle as pkl
from sklearn.linear_model import Lasso
from sklearn import linear_model
import sqlite3
from Server import Server
import sys
from DeepLearning import DeepLearning


class Regression(DeepLearning):
    def __init__(self):
        super(DeepLearning, self).__init__()
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
        self.dl = 'cnn'

    def to_server(self):
        which = 'newton'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='04:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def dl_pred(self):
        contents = {'pos': [], 'neg': []}
        for sign in contents.keys():
            fpath = os.path.join(self.root, 'Papers/Lasso', 'pairs_{}.xlsx'.format(sign))
            df = pd.read_excel(fpath, index_col=0)
            df = df[df.shape[0] // 2:]

            fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector.xlsx'.format(self.scope))
            df_gene = pd.read_excel(fpath, index_col=0)

            fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_mir_{}_score_vector.xlsx'.format(self.scope))
            df_mir = pd.read_excel(fpath, index_col=0)

            for mir in df.index:
                genes = df.loc[mir, 'GENEs'].split(';')
                for gene in genes:
                    temp = pd.concat([df_mir.loc[mir, :].round(4), df_gene.loc[gene, :].round(4)], axis=1).T
                    contents[sign].append(temp)

        X = []
        Y = []
        for sign, each in contents.items():
            if sign == 'pos':
                col = 0
            else:
                col = 1

            for i, row in enumerate(each):
                X.append(row.values)
                y = np.zeros(2)
                y[col] = 1
                Y.append(y)

        x_pred = np.dstack(X)
        x_pred = np.transpose(x_pred, (2, 0, 1))
        y_pred = np.vstack(Y)
        if self.dl == 'cnn':
            x_pred = x_pred.reshape((*x_pred.shape, 1))
        if self.dl == 'cnn':
            fname = os.path.join(self.root, 'Papers/Lasso', '{}_VGG.model'.format(self.__class__.__name__))
        elif self.dl == 'lstm':
            fname = os.path.join(self.root, 'Papers/Lasso', '{}_LSTM.model'.format(self.__class__.__name__))
        self.prediction(x_pred, y_pred, fname)

    def dl_run(self):
        contents = {'pos': [], 'neg': []}
        for sign in contents.keys():
            fpath = os.path.join(self.root, 'Papers/Lasso', 'pairs_{}.xlsx'.format(sign))
            df = pd.read_excel(fpath, index_col=0)
            df = df[: df.shape[0] // 2]

            fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector.xlsx'.format(self.scope))
            df_gene = pd.read_excel(fpath, index_col=0)

            fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_mir_{}_score_vector.xlsx'.format(self.scope))
            df_mir = pd.read_excel(fpath, index_col=0)

            for mir in df.index:
                genes = df.loc[mir, 'GENEs'].split(';')
                for gene in genes:
                    temp = pd.concat([df_mir.loc[mir, :].round(4), df_gene.loc[gene, :].round(4)], axis=1).T
                    contents[sign].append(temp)

        with open(self.__class__.__name__ + '.cha', 'wb') as f:
            pkl.dump(contents, f)

        with open(self.__class__.__name__ + '.cha', 'rb') as f:
            contents = pkl.load(f)

        X = []
        Y = []
        for sign, each in contents.items():
            if sign == 'pos':
                col = 0
            else:
                col = 1

            for i, row in enumerate(each):
                X.append(row.values)
                y = np.zeros(2)
                y[col] = 1
                Y.append(y)

        x_trn = np.dstack(X)
        x_trn = np.transpose(x_trn, (2, 0, 1))
        y_trn = np.vstack(Y)

        if self.dl == 'cnn':
            x_trn = x_trn.reshape((*x_trn.shape, 1))
        x_trn, y_trn, x_tst, y_tst = self.split_data_set(x_trn, y_trn, ratio=0.2)
        if self.dl == 'cnn':
            fname = os.path.join(self.root, 'Papers/Lasso', '{}_VGG.model'.format(self.__class__.__name__))
            self.vgg_like_convnet(x_trn, y_trn, x_tst, y_tst, fname, layer_sizes=[4, 8, 16], filter_sizes=[(1, 3), (1, 3), (1, 3)])
        elif self.dl == 'lstm':
            fname = os.path.join(self.root, 'Papers/Lasso', '{}_LSTM.model'.format(self.__class__.__name__))
            self.sequence_LSTM(x_trn, y_trn, x_tst, y_tst, fname)

    # def run(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_corr2.xlsx'.format(self.scope))
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector.xlsx'.format(self.scope))
    #     df_gene = pd.read_excel(fpath, index_col=0)
    #
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_mir_{}_score_vector.xlsx'.format(self.scope))
    #     df_mir = pd.read_excel(fpath, index_col=0)
    #
    #     contents = {}
    #     for mir in df.index:
    #         genes = df.loc[mir, 'GENEs'].split(';')
    #         corr = df.loc[mir, 'corr (pearson)']
    #
    #         # mir_values = ';'.join(df_mir.loc[mir, :].round(4).astype(str))
    #         for gene in genes:
    #             if corr not in contents:
    #                 # corr = corr + ':' + df.loc[mir, 'GENEs'] + ':' + mir
    #                 corr = ';'.join(df_mir.loc[mir, :].round(4).astype(str)) + ':' + df.loc[mir, 'GENEs'] + ':' + mir
    #                 contents[corr] = [df_gene.loc[gene, :].values]
    #             else:
    #                 contents[corr].append(df_gene.loc[gene, :].values)
    #
    #     with open(self.__class__.__name__ + '.cha', 'wb') as f:
    #         pkl.dump(contents, f)
    #
    #     with open(self.__class__.__name__ + '.cha', 'rb') as f:
    #         contents = pkl.load(f)
    #
    #     data = []
    #     columns = ['miRNA', 'gene', 'coeff', 'intercept']
    #     clf = linear_model.Lasso(alpha=0.1)
    #     for corr, genes in contents.items():
    #         if len(genes) < 2:
    #             continue
    #         corr, gene_names, mir = np.array(corr.split(':'))
    #         corr = np.array(corr.split(';')).astype(float)
    #         gene_names = gene_names.split(';')
    #
    #         genes = np.stack(genes).T
    #         clf.fit(genes, corr)
    #         Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
    #               normalize=False, positive=False, precompute=False, random_state=None,
    #               selection='cyclic', tol=1e-4, warm_start=False)
    #         for c, gname in zip(clf.coef_, gene_names):
    #             data.append([mir, gname, c, clf.intercept_])
    #
    #     df = pd.DataFrame(data, columns=columns)
    #     out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', self.__class__.__name__ + '.xlsx')
    #     df.to_excel(out_path, index=None)

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_corr2.xlsx'.format(self.scope))
        df = pd.read_excel(fpath, index_col=0)

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_{}_score_vector.xlsx'.format(self.scope))
        df_gene = pd.read_excel(fpath, index_col=0)

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'human_cell_line_hCAGE_mir_{}_score_vector.xlsx'.format(self.scope))
        df_mir = pd.read_excel(fpath, index_col=0)

        contents = {}
        for mir in df.index:
            genes = df.loc[mir, 'GENEs'].split(';')
            corr = df.loc[mir, 'corr']

            # mir_values = ';'.join(df_mir.loc[mir, :].round(4).astype(str))
            for gene in genes:
                if corr not in contents:
                    # corr = corr + ':' + df.loc[mir, 'GENEs'] + ':' + mir
                    corr = ';'.join(df_mir.loc[mir, :].round(4).astype(str)) + ':' + df.loc[mir, 'GENEs'] + ':' + mir
                    contents[corr] = [df_gene.loc[gene, :].values]
                else:
                    contents[corr].append(df_gene.loc[gene, :].values)

        with open(self.__class__.__name__ + '.cha', 'wb') as f:
            pkl.dump(contents, f)

        with open(self.__class__.__name__ + '.cha', 'rb') as f:
            contents = pkl.load(f)

        data = []
        columns = ['miRNA', 'gene', 'coeff', 'intercept']
        clf = linear_model.Lasso(alpha=0.1)
        for corr, genes in contents.items():
            corr, gene_names, mir = np.array(corr.split(':'))
            corr = np.array(corr.split(';')).astype(float)
            gene_names = gene_names.split(';')
            for gname, gene in zip(gene_names, genes):
                gene = gene.reshape((len(gene), 1))
                clf.fit(gene, corr)
                Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
                   normalize=False, positive=False, precompute=False, random_state=None,
                   selection='cyclic', tol=1e-4, warm_start=False)
                data.append([mir, gname, *clf.coef_, clf.intercept_])

        df = pd.DataFrame(data, columns=columns)
        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', self.__class__.__name__ + '.xlsx')
        df.to_excel(out_path, index=None)

    def histogram(self):
        import matplotlib.pyplot as plt

        df = pd.read_excel(self.__class__.__name__ + '.xlsx')
        hist, bin = np.histogram(df['coeff'].values, bins=100)
        plt.bar(bin[:-1], hist)
        plt.show()

    def filter_by_lasso(self):
        fname = self.__class__.__name__
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines', fname + '.xlsx')
        df = pd.read_excel(fpath)
        df_grp = df.groupby('miRNA')
        print(len(df_grp))

        contents = []
        for mir, df_sub in df_grp:
            df_sub = df_sub.sort_values(by=['coeff'], ascending=False)
            df_sub = df_sub[df_sub['coeff'].abs() > 0.5]
            if df_sub.empty:
                continue
            contents.append([mir, ';'.join(df_sub['gene'], ), ';'.join(df_sub['coeff'].round(2).astype(str))])

        df_res = pd.DataFrame(data=contents, columns=['miRNA', 'GENEs', 'Coeffs'])
        out_path = os.path.join(self.root, 'database/Fantom/v5/cell_lines', fname + '_filter.xlsx')
        df_res.to_excel(out_path, index=None)

    def get_distance(self):
        fpath = self.__class__.__name__ + '.xlsx'
        df = pd.read_excel(self.__class__.__name__ + '.xlsx')

        mir_path = os.path.join(self.root, 'database', 'fantom5.db')
        con_mir = sqlite3.connect(mir_path)

        gene_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
        con_gene = sqlite3.connect(gene_path)

        appendex = []
        for idx in df.index:
            if idx % 1000 == 0 or idx + 1 == df.shape[0]:
                print('{:0.2f}%'.format(100 * (idx + 1) / df.shape[0]))

            mirna = df.loc[idx, 'miRNA']
            gene = df.loc[idx, 'gene']

            df_mir = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE premiRNA='{}'".format(mirna), con_mir, index_col='premiRNA')
            df_gene = pd.read_sql_query("SELECT * FROM 'human_genes_wo_duplicates' WHERE gene_name='{}'".format(gene), con_gene, index_col='gene_name')
            if df_mir.empty or df_gene.empty:
                appendex.append([None, None])
                continue

            mir_chr = df_mir.loc[mirna, 'chromosome']
            mir_strand = df_mir.loc[mirna, 'strand']
            gene_chr = df_gene.loc[gene, 'chromosome']
            gene_strand = df_gene.loc[gene, 'strand']

            if mir_chr == gene_chr and mir_strand == gene_strand:
                if mir_strand == '+':
                    offset = 'start'
                else:
                    offset = 'end'
                distance = abs(df_gene.loc[gene, offset] - df_mir.loc[mirna, offset])
                appendex.append([mir_chr, distance])
            else:
                appendex.append([mir_chr, None])
        df_app = pd.DataFrame(appendex)
        df = pd.concat([df, df_app], axis=1)
        df.to_excel(fpath, index=None)

    def add_truth(self):
        fpath_ref = os.path.join(self.root, 'database', 'miRTarBase.db')
        con_ref = sqlite3.connect(fpath_ref)
        df_ref = pd.read_sql_query("SELECT * FROM 'miRTarBase' WHERE Species_miRNA='Homo sapiens'", con_ref)
        df_ref['miRNA'] = df_ref['miRNA'].str.lower()
        df_ref = df_ref.set_index('miRNA')

        fpath_fan = os.path.join(self.root, 'database', 'fantom5.db')
        con_fan = sqlite3.connect(fpath_fan)
        df_fan = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates'", con_fan, index_col='premiRNA')

        df_cor = pd.read_excel('Regression.xlsx', index_col=0)
        df_cor.loc[:, 'index'] = np.arange(df_cor.shape[0])
        df_cor.loc[:, 'truth'] = False

        dfs = []
        for pmir in set(df_cor.index):
            df_cor_sub = df_cor.loc[pmir]
            mir = df_fan.loc[pmir, 'miRNA'].lower()
            gene = df_cor_sub['gene']
            if mir not in df_ref.index:
                continue

            tgenes = df_ref.loc[mir, 'Target Gene']
            common = list(set.intersection(set(gene), set(tgenes)))
            if len(common) > 0:
                df_cor_sub = df_cor_sub.set_index('gene', drop=False)
                idx = df_cor_sub.loc[common, 'index'].values
                df_cor_sub = df_cor_sub.set_index('index', drop=False)
                df_cor_sub.loc[idx, 'truth'] = True
                df_cor_sub.loc[:, 'miRNA'] = pmir
                dfs.append(df_cor_sub)

        df_cor = pd.concat(dfs)
        df_cor.to_excel(self.__class__.__name__ + '2.xlsx')

    def comparison(self):
        df = pd.read_excel(self.__class__.__name__ + '.xlsx', index_col=0)

        fpath = os.path.join(self.root, 'database', 'miRTarBase.db')
        con = sqlite3.connect(fpath)

        mir_path = os.path.join(self.root, 'database', 'fantom5.db')
        con_mir = sqlite3.connect(mir_path)

        dfs = []
        df_ref = pd.read_sql_query("SELECT * FROM 'miRTarBase' WHERE Species_miRNA='Homo sapiens'", con)
        for i, pmir in enumerate(df.index):
            if i % 1000 == 0 or i + 1 == df.shape[0]:
                print('{:0.2f}%'.format(100 * (i + 1) / df.shape[0]))

            df_sub = df.loc[pmir, :]

            mir = pd.read_sql_query("SELECT miRNA FROM 'human_promoters_wo_duplicates' WHERE premiRNA='{}'".format(pmir), con_mir)
            if len(mir) == 0:
                continue

            df_ref_sub = df_ref[df_ref['miRNA'] == mir.loc[0, 'miRNA']]

            df_sub = df_sub.reset_index()
            df_sub = df_sub.set_index('gene')
            gene_inter = list(set.intersection(set(df_sub.index), set(df_ref_sub['Target Gene'])))
            if len(gene_inter) == 0:
                continue

            df_sub = df_sub.loc[gene_inter, :]
            dfs.append(df_sub.reset_index())

        df_res = pd.concat(dfs)
        df_res.to_excel(self.__class__.__name__ + '2.xlsx', index=None)

    #
    # def run(self):
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'human_cell_line_hCAGE_{}_score_vector_corr2.xlsx'.format(self.scope))
    #     df = pd.read_excel(fpath, index_col=0)
    #
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'human_cell_line_hCAGE_{}_score_vector.xlsx'.format(self.scope))
    #     df_gene = pd.read_excel(fpath, index_col=0)
    #
    #     fpath = os.path.join(self.root, 'database/Fantom/v5/tissues/out', 'human_cell_line_hCAGE_mir_{}_score_vector.xlsx'.format(self.scope))
    #     df_mir = pd.read_excel(fpath, index_col=0)
    #
    #     contents = {}
    #     for mir in df.index:
    #         genes = df.loc[mir, 'GENEs'].split(';')
    #         # corr = df.loc[mir, 'corr (pearson)']
    #
    #         mir_values = ';'.join(df_mir.loc[mir, :].round(4).astype(str))
    #         for gene in genes:
    #             if mir_values in contents:
    #                 contents[mir_values].append(df_gene.loc[gene, :].round(4).values)
    #             else:
    #                 contents[mir_values] = [df_gene.loc[gene, :].round(4).values]
    #
    #     with open(self.__class__.__name__ + '.cha', 'wb') as f:
    #         pkl.dump(contents, f)
    #
    #     with open(self.__class__.__name__ + '.cha', 'rb') as f:
    #         contents = pkl.load(f)
    #
    #     data = []
    #     sheet = []
    #     clf = linear_model.Lasso(alpha=0.1)
    #     for mir_values, gene_values in contents.items():
    #         mir_values = np.array(mir_values.split(';')).astype(float).round(4)
    #         gene_values = np.array(gene_values).T
    #         clf.fit(gene_values, mir_values)
    #
    #         Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
    #            normalize=False, positive=False, precompute=False, random_state=None,
    #            selection='cyclic', tol=1e-4, warm_start=False)
    #         data.append([*clf.coef_, clf.intercept_])
    #         sheet.append(clf.intercept_)
    #
    #     df = pd.Series(sheet, index=df.index)
    #     df.to_excel(self.__class__.__name__ + '.xlsx')
    #     with open(self.__class__.__name__ + '.cha', 'wb') as f:
    #         pkl.dump(data, f)

    def filter_by_db(self):
        fname = 'human_cell_line_hCAGE_{}_score_corr2_with_target_indicator_ts.csv'.format(self.scope)
        df_ref = pd.read_csv(os.path.join(self.root, 'Papers/Lasso', fname), index_col=0)

        pairs = {'pos': {}, 'neg': {}}
        for mir in df_ref.index:
            genes = df_ref.loc[mir, 'GENEs'].split(';')
            ti = map(int, df_ref.loc[mir, 'Target Indicator'].split(';'))
            for g, i in zip(genes, ti):
                if i > 0:
                    if mir not in pairs['pos']:
                        pairs['pos'][mir] = [g]
                    else:
                        pairs['pos'][mir].append(g)
                else:
                    if mir not in pairs['neg']:
                        pairs['neg'][mir] = [g]
                    else:
                        pairs['neg'][mir].append(g)
        return pairs

    def evaluation(self):
        pairs = self.filter_by_db()

        for sign, pair in pairs.items():
            contents = []
            for mir, genes in pair.items():
                contents.append([mir, ';'.join(genes)])
            df = pd.DataFrame(data=contents, columns=['miRNA', 'GENEs'])
            out_path = os.path.join(self.root, 'Papers/Lasso', 'pairs_{}.xlsx'.format(sign))
            df.to_excel(out_path, index=None)

        # pidx = []
        # for idx in df_reg.index:
        #     mir = df_reg.loc[idx, 'miRNA']
        #     gene = df_reg.loc[idx, 'gene']
        #     if mir in pairs and gene in pairs[mir]:
        #         pidx.append(idx)
        #
        # nidx = list(set(df_reg.index) - set(pidx))
        # df_reg_pos = df_reg.loc[pidx]
        # df_reg_neg = df_reg.loc[nidx]
        #
        # df_reg_pos.to_excel(fpath.replace('.xlsx', '_pos.xlsx'), index=None)
        # df_reg_neg.to_excel(fpath.replace('.xlsx', '_neg.xlsx'), index=None)


if __name__ == '__main__':
    rg = Regression()
    if rg.hostname == 'mingyu-Precision-Tower-7810':
        # rg.to_server()
        rg.run()
        rg.filter_by_lasso()
        # rg.dl_pred()
        # rg.evaluation()
    else:
        rg.run()
