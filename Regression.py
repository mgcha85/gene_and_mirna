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
from Database import Database
import matplotlib.pyplot as plt


class Regression(DeepLearning):
    def __init__(self, root):
        super(DeepLearning, self).__init__()
        self.root = root
        self.table_names = {}

        self.scope = 500
        self.dl = 'cnn'

    def to_server(self):
        from Server import Server
        import sys

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
            corr = df.loc[mir, 'corr (pearson)']

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

    # def get_distance(self):
    #     fpath = self.__class__.__name__ + '.xlsx'
    #     df = pd.read_excel(self.__class__.__name__ + '.xlsx')
    #
    #     mir_path = os.path.join(self.root, 'database', 'fantom5.db')
    #     con_mir = sqlite3.connect(mir_path)
    #
    #     gene_path = os.path.join(self.root, 'database/gencode', 'gencode.v30lift37.basic.annotation.db')
    #     con_gene = sqlite3.connect(gene_path)
    #
    #     appendex = []
    #     for idx in df.index:
    #         if idx % 1000 == 0 or idx + 1 == df.shape[0]:
    #             print('{:0.2f}%'.format(100 * (idx + 1) / df.shape[0]))
    #
    #         mirna = df.loc[idx, 'miRNA']
    #         gene = df.loc[idx, 'gene']
    #
    #         df_mir = pd.read_sql("SELECT * FROM 'human_promoters_wo_duplicates' WHERE premiRNA='{}'".format(mirna), con_mir, index_col='premiRNA')
    #         df_gene = pd.read_sql("SELECT * FROM 'human_genes_wo_duplicates' WHERE gene_name='{}'".format(gene), con_gene, index_col='gene_name')
    #         if df_mir.empty or df_gene.empty:
    #             appendex.append([None, None])
    #             continue
    #
    #         mir_chr = df_mir.loc[mirna, 'chromosome']
    #         mir_strand = df_mir.loc[mirna, 'strand']
    #         gene_chr = df_gene.loc[gene, 'chromosome']
    #         gene_strand = df_gene.loc[gene, 'strand']
    #
    #         if mir_chr == gene_chr and mir_strand == gene_strand:
    #             if mir_strand == '+':
    #                 offset = 'start'
    #             else:
    #                 offset = 'end'
    #             distance = abs(df_gene.loc[gene, offset] - df_mir.loc[mirna, offset])
    #             appendex.append([mir_chr, distance])
    #         else:
    #             appendex.append([mir_chr, None])
    #     df_app = pd.DataFrame(appendex)
    #     df = pd.concat([df, df_app], axis=1)
    #     df.to_excel(fpath, index=None)

    def add_truth(self):
        fpath_ref = os.path.join(self.root, 'database', 'miRTarBase.db')
        con_ref = sqlite3.connect(fpath_ref)
        df_ref = pd.read_sql("SELECT * FROM 'miRTarBase' WHERE Species_miRNA='Homo sapiens'", con_ref)
        df_ref['miRNA'] = df_ref['miRNA'].str.lower()
        df_ref = df_ref.set_index('miRNA')

        fpath_fan = os.path.join(self.root, 'database', 'fantom5.db')
        con_fan = sqlite3.connect(fpath_fan)
        df_fan = pd.read_sql("SELECT * FROM 'human_promoters_wo_duplicates'", con_fan, index_col='premiRNA')

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
        df_ref = pd.read_sql("SELECT * FROM 'miRTarBase' WHERE Species_miRNA='Homo sapiens'", con)
        for i, pmir in enumerate(df.index):
            if i % 1000 == 0 or i + 1 == df.shape[0]:
                print('{:0.2f}%'.format(100 * (i + 1) / df.shape[0]))

            df_sub = df.loc[pmir, :]

            mir = pd.read_sql("SELECT miRNA FROM 'human_promoters_wo_duplicates' WHERE premiRNA='{}'".format(pmir), con_mir)
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

    def merge_table(self, fpath):
        from Database import Database

        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        dfs = []
        for tname in tlist:
            dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
        df = pd.concat(dfs)
        df = df.set_index(df.columns[0])
        df.index.name = df.columns[0]
        return df

    def get_trn_data(self, hbw, type):
        fpaths = {'mir': os.path.join(self.root, 'database/Fantom/v5/{}'.format(type), 'sum_fan_mir_{}.db'.format(hbw)),
                  'gene': os.path.join(self.root, 'database/Fantom/v5/{}'.format(type), 'sum_fan_gene_{}.db'.format(hbw))}
        dfs = {}
        for label, fpath in fpaths.items():
            dfs[label] = self.merge_table(fpath)
            print('[{}] #:{}'.format(label, dfs[label].shape[0]))

        if type == 'cell_lines':
            gene = dfs['gene'].loc[:, '10964C':].T
            mir = dfs['mir'].loc[:, '10964C':].T
        elif type == 'tissues':
            gene = dfs['gene'].loc[:, 'appendix':].T
            mir = dfs['mir'].loc[:, 'appendix':].T
        else:
            return
        return gene, mir

    def regression(self, hbw, opt, type='cell_lines'):
        gene, mir = self.get_trn_data(hbw, type)

        # # miRNAs in RNA-seq
        # fpath = os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_mir.db')
        # con = sqlite3.connect(fpath)
        # df_rna = pd.read_sql("SELECT * FROM 'MIR_expression'", con, index_col='Name')
        # mir = mir[df_rna.index]

        gene = gene - gene.mean(axis=0)
        mir = mir - mir.mean(axis=0)

        # gene = (gene - gene.mean(axis=0)) / gene.std(axis=0)
        # gene = gene.dropna(how='any', axis=1)
        # mir = (mir - mir.mean(axis=0)) / mir.std(axis=0)
        # mir = mir.dropna(how='any', axis=1)

        clf = linear_model.Lasso(alpha=0.1)
        clf.fit(gene, mir)
        Lasso(alpha=0.1, copy_X=True, fit_intercept=False, max_iter=1000,
              normalize=False, positive=False, precompute=False, random_state=None,
              selection='cyclic', tol=1e-3, warm_start=False)

        df_coef = pd.DataFrame(clf.coef_, index=mir.columns, columns=gene.columns).T
        df_inter = pd.Series(clf.intercept_, index=mir.columns)

        fpath = os.path.join(self.root, 'database/Fantom/v5/{}/out'.format(type), 'regression_{}_{}.db'.format(hbw, opt))
        con = sqlite3.connect(fpath)

        gene = gene.T
        gene.index.name = 'tid'
        gene.to_sql('X', con, if_exists='replace')

        mir = mir.T
        mir.index.name = 'miRNA'
        mir.to_sql('Y', con, if_exists='replace')
        df_coef.index.name = 'transcript_id'
        df_coef.to_sql('coefficient', con, if_exists='replace')
        df_inter.to_sql('intercept', con, if_exists='replace')
        print('[{}] {:0.4f}'.format(hbw, clf.score(gene.T, mir.T)))

    def regression_rna(self):
        fpaths = {'mir': os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_mir.db'),
                  'gene': os.path.join(self.root, 'database/RNA-seq/out', 'RNA_seq_gene.db')}
        dfs = {}
        for label, fpath in fpaths.items():
            df = pd.read_sql("SELECT * FROM '{}_expression'".format(label), sqlite3.connect(fpath))
            dfs[label] = df.set_index(df.columns[0], drop=True)
            print('[{}] #:{}'.format(label, dfs[label].shape[0]))

        clf = linear_model.Lasso(alpha=0.1)
        gene = dfs['gene'].T
        mir = dfs['mir'].T

        gene = gene - gene.mean(axis=0)
        mir = mir - mir.mean(axis=0)

        # gene = (gene - gene.mean(axis=0)) / gene.std(axis=0)
        # gene = gene.dropna(how='any', axis=1)
        # mir = (mir - mir.mean(axis=0)) / mir.std(axis=0)
        # mir = mir.dropna(how='any', axis=1)

        clf.fit(gene, mir)
        Lasso(alpha=0.1, copy_X=True, fit_intercept=False, max_iter=1000,
              normalize=False, positive=False, precompute=False, random_state=None,
              selection='cyclic', tol=1e-3, warm_start=False)

        df_coef = pd.DataFrame(clf.coef_, index=mir.columns, columns=gene.columns).T
        df_inter = pd.Series(clf.intercept_, index=mir.columns)

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_rna.db')
        con = sqlite3.connect(fpath)

        gene.T.to_sql('X', con, if_exists='replace')
        mir.T.to_sql('Y', con, if_exists='replace')
        df_coef.to_sql('coefficient', con, if_exists='replace')
        df_inter.to_sql('intercept', con, if_exists='replace')

    def compare(self):
        from scipy.stats import hypergeom
        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100.db')
        con_fan = sqlite3.connect(fpath_fan)

        fpath_rna = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_rna.db')
        con_rna = sqlite3.connect(fpath_rna)

        df_fan = pd.read_sql("SELECT * FROM 'result'", con_fan, index_col='miRNA')
        df_rna = pd.read_sql("SELECT * FROM 'result'", con_rna, index_col='miRNA')

        mirList = set.intersection(set(df_fan.index), set(df_rna.index))
        contents = []
        for mir in mirList:
            gene_fan = set(df_fan.loc[mir, 'gene_name'].split(';'))
            gene_rna = set(df_rna.loc[mir, 'gene_name'].split(';'))
            intersection = set.intersection(gene_fan, gene_rna)

            U = set.union(gene_fan, gene_rna)
            m = len(U - gene_fan)
            n = len(gene_fan)

            k = len(gene_rna)
            q = len(intersection)
            p = hypergeom.sf(q-1, n+m, m, k)

            contents.append([mir, ';'.join(sorted(list(gene_fan))), ';'.join(sorted(list(gene_rna))), m, n, q, p])
        df_res = pd.DataFrame(contents, columns=['miRNA', 'genes (fantom)', 'genes (RNA-seq)', 'm', 'n', 'i', 'p-value'])
        df_res.to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_compare.xlsx'), index=None)

    def cross_stats(self, hbw, opt, type='cell_lines'):
        dirname = os.path.join(self.root, 'database/Fantom/v5/{}/out/cross_validation'.format(type), opt)
        flist = os.listdir(dirname)
        flist = sorted([x for x in flist if x.endswith('.db') and '_test' in x])
        columns__ = ['miRNA', 'median_diff', 'mean_diff', 'std_diff', 'median_expr', 'mean_expr', 'std_expr', 'med_diff_expr_ratio']

        writer = pd.ExcelWriter(os.path.join(dirname, 'cross_stats.xlsx'), engine='xlsxwriter')
        summary, sig_mir, num_sig = [], [], []

        for fname in flist:
            report_trn, report_test = [], []
            # test
            fpath_test = os.path.join(dirname, fname)
            con_test = sqlite3.connect(fpath_test)
            # train
            fpath_trn = os.path.join(dirname, fname.replace('_test', '_train'))
            con_trn = sqlite3.connect(fpath_trn)

            # load test data
            # B_test = pd.read_sql("SELECT * FROM 'coefficient'", con_test, index_col='tid')
            X_test = pd.read_sql("SELECT * FROM 'X'", con_test, index_col='tid')
            Y_test = pd.read_sql("SELECT * FROM 'Y'", con_test, index_col='miRNA')

            # load train data
            B_trn = pd.read_sql("SELECT * FROM 'coefficient'", con_trn, index_col='tid').T
            X_trn = pd.read_sql("SELECT * FROM 'X'", con_trn, index_col='tid')
            Y_trn = pd.read_sql("SELECT * FROM 'Y'", con_trn, index_col='miRNA')

            # matrix multiplication
            Yh_trn = np.dot(B_trn, X_trn)
            E_trn = (Y_trn - Yh_trn).abs()

            # matrix multiplication with plugin
            Yh_test = np.dot(B_trn, X_test)
            E_test = (Y_test - Yh_test).abs()

            # get statistics
            df_rep = []
            for E, report, Y, label in zip([E_test, E_trn], [report_test, report_trn], [Y_test, Y_trn], [' (test)', ' (train)']):
                for idx in E.index:
                    med_diff = E.loc[idx].median()
                    mean_diff = E.loc[idx].mean()
                    std_diff = E.loc[idx].std()

                    med_expr = Y.loc[idx].median()
                    mean_expr = Y.loc[idx].mean()
                    std_expr = Y.loc[idx].std()
                    if med_expr == 0:
                        med_diff_expr_ratio = float('inf')
                    else:
                        med_diff_expr_ratio = med_diff / med_expr
                    report.append([idx, med_diff, mean_diff, std_diff, med_expr, mean_expr, std_expr, med_diff_expr_ratio])

                # plt.hist(E.values.flatten(), bins='auto')
                # plt.xlabel('Error')
                # plt.ylabel('Frequency')
                # plt.title(os.path.splitext(fname)[0])
                # plt.savefig(os.path.join(dirname, os.path.splitext(fname)[0] + '.png'))
                # plt.close()

                columns = []
                for i in range(len(columns__)):
                    columns.append(columns__[i] + label)
                df_rep.append(pd.DataFrame(report, columns=columns))

            df_rep = pd.concat(df_rep, axis=1)
            df_rep['med_der_ratio'] = df_rep['med_diff_expr_ratio (test)'] / df_rep['med_diff_expr_ratio (train)']
            df_rep['median_diff_ratio'] = (df_rep['median_diff (test)'].abs() - df_rep['median_diff (train)']).abs()
            # print(len(df_rep[df_rep['med_der_ratio'] < 10]))
            num_sig.append(len(df_rep[df_rep['median_diff_ratio'] < 1.1]))

            _, hbw, num, type = os.path.splitext(fname)[0].split('_')
            sig_mir.append(set(df_rep[df_rep['median_diff_ratio'] < 1.1]['miRNA (test)']))
            df_rep.sort_values(by='median_diff_ratio').to_excel(writer, sheet_name=num, index=None)
            summary.append(df_rep.mean())

        sig_mir = set.intersection(*sig_mir)
        with open('singificant_miRNA.txt', 'wt') as f:
            f.write('\n'.join(sorted(list(sig_mir))))

        with open('num_sig.txt', 'wt') as f:
            f.write('\n'.join(map(str, num_sig)))

        # write report
        df_summary = pd.concat(summary, axis=1).T
        df_summary.index = flist
        df_summary.to_excel(os.path.join(dirname, 'cross_stats_summary.xlsx'))
        writer.save()
        writer.close()

    def compare_tissue_cross(self, opt):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', opt)
        fpath = os.path.join(dirname, 'cross_stats.xlsx')
        reader = pd.ExcelFile(fpath)
        snames = ['tissue'] + reader.sheet_names

        df_tis = pd.read_excel(os.path.join(self.root, 'database/Fantom/v5/tissues/out/cross_validation', opt, 'cross_stats.xlsx'))
        df_res = pd.DataFrame(index=snames, columns=snames)

        dfs = {'tissue': df_tis}
        for i in snames[1:]:
            dfs[i] = reader.parse(i)

        for s1 in snames:
            for s2 in snames:
                if s1 == s2:
                    continue
                df1 = dfs[s1]
                df2 = dfs[s2]
                mir1 = set(df1[df1['median_diff_ratio'] < 1.1]['miRNA (test)'])
                mir2 = set(df2[df2['median_diff_ratio'] < 1.1]['miRNA (test)'])
                mir = set.intersection(mir1, mir2)
                df_res.loc[s1, s2] = len(mir) / len(set.union(mir1, mir2))
        df_res.to_excel(os.path.join(dirname, 'cross_summary.xlsx'))

    def check_overlap(self, opt, type):
        dirname = os.path.join(self.root, 'database/Fantom/v5/{}/out/cross_validation'.format(type), opt)
        fpath = os.path.join(dirname, 'cross_stats.xlsx')
        reader = pd.ExcelFile(fpath)
        snames = reader.sheet_names

        df_res = pd.DataFrame(index=snames, columns=snames)
        for i in snames:
            df1 = reader.parse(i)
            for j in snames:
                if i == j:
                    continue
                df2 = reader.parse(j)
                mir1 = set(df1[df1['median_diff_ratio'] < 1.1]['miRNA (test)'])
                mir2 = set(df2[df2['median_diff_ratio'] < 1.1]['miRNA (test)'])
                mir = set.intersection(mir1, mir2)
                df_res.loc[i, j] = len(mir) / df1.shape[0]
        df_res.to_excel(os.path.join(dirname, 'cross_summary.xlsx'))

    def cross_regression(self, hbw, opt, N=10, type='cell_lines'):
        fpaths = {'mir': os.path.join(self.root, 'database/Fantom/v5/{}'.format(type), 'sum_fan_mir_{}.db'.format(hbw)),
                  'gene': os.path.join(self.root, 'database/Fantom/v5/{}'.format(type), 'sum_fan_gene_{}.db'.format(hbw))}
        dfs = {}
        for label, fpath in fpaths.items():
            dfs[label] = self.merge_table(fpath)
            print('[{}] #:{}'.format(label, dfs[label].shape[0]))

        clf = linear_model.Lasso(alpha=0.1)
        if type == 'cell_lines':
            df_gene = dfs['gene'].loc[:, '10964C':].T
            df_mir = dfs['mir'].loc[:, '10964C':].T
        elif type == 'tissues':
            df_gene = dfs['gene'].loc[:, 'achilles tendon':].T
            df_mir = dfs['mir'].loc[:, 'achilles tendon':].T
        else:
            return

        def get_batch(df, i, axis=0):
            if axis == 0:
                n = df.shape[0]
            elif axis == 1:
                n = df.shape[1]
            else:
                raise Exception("wrong axis")

            m = (n + N - 1) // N
            idx = np.arange(m * i, m * (i + 1))
            idx = idx[idx < n]
            dixd = list(df.index[idx])
            if axis == 0:
                trn, test = df.drop(dixd), df.loc[dixd, :]
            elif axis == 1:
                trn, test = df.drop(dixd, axis=1), df.loc[:, dixd]
            return trn, test

        for i in range(N-1, N):
            gene_trn, gene_test = get_batch(df_gene, i)
            mir_trn, mir_test = get_batch(df_mir, i)

            for gt, mt, label in zip([gene_trn, gene_test], [mir_trn, mir_test], ['train', 'test']):
                # normalization
                gt -= gt.mean(axis=0)
                mt -= mt.mean(axis=0)
                # gt = (gt - gt.std(axis=0)) / gt.mean(axis=0)
                # mt = (mt - mt.std(axis=0)) / mt.mean(axis=0)

                clf.fit(gt, mt)
                Lasso(alpha=0.1, copy_X=True, fit_intercept=False, max_iter=1000,
                      normalize=False, positive=False, precompute=False, random_state=None,
                      selection='cyclic', tol=1e-3, warm_start=False)

                df_coef = pd.DataFrame(clf.coef_, index=mt.columns, columns=gt.columns).T
                df_inter = pd.Series(clf.intercept_, index=mt.columns)

                dirname = os.path.join(self.root, 'database/Fantom/v5/{}/out/cross_validation'.format(type), opt)
                if not os.path.exists(dirname):
                    os.mkdir(dirname)
                fpath = os.path.join(dirname, 'regression_{}_{}_{}.db'.format(hbw, i, label))
                con = sqlite3.connect(fpath)

                gt = gt.T
                gt.index.name = 'tid'
                gt.to_sql('X', con, if_exists='replace')

                mt = mt.T
                mt.index.name = 'miRNA'
                mt.to_sql('Y', con, if_exists='replace')
                df_coef.index.name = 'tid'
                df_coef.to_sql('coefficient', con, if_exists='replace')
                df_inter.index.name = 'miRNA'
                df_inter.to_sql('intercept', con, if_exists='replace')
                print('[{}] {:0.4f}'.format(hbw, clf.score(gt.T, mt.T)))

    def get_plugin_distance(self, opt):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', opt)

        contents = {'test': [], 'train': []}
        index = []
        for fname in sorted([x for x in os.listdir(dirname) if x.endswith('.db')]):
            if 'train' in fname:
                contents['train'].append(fname)
            elif 'test' in fname:
                contents['test'].append(fname)
            else:
                continue

        report = []
        for fname in sorted(contents['train']):
            name, tpm, num, type = os.path.splitext(fname)[0].split('_')
            index.append(num)
            con_trn = sqlite3.connect(os.path.join(dirname, fname))
            fname = '_'.join([name, tpm, num, 'test']) + '.db'
            con_tst = sqlite3.connect(os.path.join(dirname, fname))

            B = pd.read_sql("SELECT * FROM 'coefficient'", con_trn, index_col='tid')
            X = pd.read_sql("SELECT * FROM 'X'", con_tst, index_col='tid')
            Y = pd.read_sql("SELECT * FROM 'Y'", con_tst, index_col='miRNA')

            Yh = np.dot(B.T, X)
            distance = (Y - Yh).abs().values.sum() / Y.size
            report.append(distance)
        print(report)
        pd.DataFrame(report, index=index).to_excel(os.path.join(dirname, 'plugin_distance.xlsx'))

    def get_distance(self, opt):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/cross_validation', opt)

        contents = {'test': [], 'train': []}
        index = []
        for fname in sorted([x for x in os.listdir(dirname) if x.endswith('.db')]):
            name, tpm, num, type = os.path.splitext(fname)[0].split('_')
            fpath = os.path.join(dirname, fname)
            con = sqlite3.connect(fpath)
            y = pd.read_sql("SELECT * FROM 'Y'", con, index_col='miRNA')
            x = pd.read_sql("SELECT * FROM 'X'", con, index_col='tid')
            coeff = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='tid')
            yh = np.dot(x, coeff)
            distance = (y - yh).abs().values.sum() / y.size
            print('y.size: {}, distance: {}, type: {}'.format(y.size, distance, type))
            contents[type].append(distance)
            if type == 'test':
                index.append('_'.join([type, num]))
        print(contents)
        pd.DataFrame(contents, index=index).to_excel(os.path.join(dirname, 'distances.xlsx'))

    def report(self, hbw, opt='nz', type='cell_lines'):
        fpath = os.path.join(self.root, 'database/Fantom/v5/{}/out'.format(type), 'regression_{}_{}.db'.format(hbw, opt))
        con = sqlite3.connect(fpath)
        # df = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='transcript_name')
        df = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='transcript_id')

        contents = []
        for col in df.columns:
            ser = df[col]
            if opt == 'nz':
                ser = ser[ser != 0]
            elif opt == 'neg':
                ser = ser[ser < 0]
            elif opt == 'pos':
                ser = ser[ser > 0]
            contents.append([col, ';'.join(ser.index)])
        pd.DataFrame(data=contents, columns=['miRNA', 'Transcripts']).to_sql('result', con, if_exists='replace', index=False)

    def add_gene_name(self, hbw, opt, type='cell_lines'):
        fpath = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        dfs = []
        for tname in tlist:
            dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='transcript_id'))
        df = pd.concat(dfs)

        res_con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/{}/out'.format(type), 'regression_{}_{}.db'.format(hbw, opt)))
        df_res = pd.read_sql("SELECT * FROM 'result' WHERE Transcripts>''", res_con)
        df_res['gene_name'] = None

        for idx in df_res.index:
            if idx % 100 == 0 or idx + 1 == df_res.shape[0]:
                print('{:0.2f}%'.format(100 * (idx + 1) / df_res.shape[0]))
            tr = df_res.loc[idx, 'Transcripts']
            if not isinstance(tr, str):
                continue

            tid = []
            gnames = []
            for t in tr.split(';'):
                if t in df.index:
                    gnames.append(df.loc[t, 'gene_name'])
                    tid.append(t)

            df_res.loc[idx, 'Transcripts'] = ';'.join(tid)
            df_res.loc[idx, 'gene_name'] = ';'.join(gnames)
        df_res.to_sql('result', res_con, if_exists='replace', index=False)

    def move(self):
        out_con = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'predictions_processed_result.db'))
        res_con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_rna.db'))
        df = pd.read_sql("SELECT * FROM 'result'", res_con, index_col='miRNA')
        df.index.name = 'pre-miRNA'
        df = df.rename(columns={'gene_name': 'genes'})
        df = df.drop('Transcripts', axis=1)
        df.to_sql('rnaseq_pre', out_con, if_exists='replace')

    def filtering(self, hbw):
        fpath = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_{}.db'.format(hbw))
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        dfs = []
        for tname in tlist:
            dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='transcript_id'))
            # dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='transcript_name'))
        df_ref = pd.concat(dfs)
        gene_names = set(df_ref['gene_name'])
        print(len(gene_names))

        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)

        tnames = ['miRTartBase_hsa', 'target_scan_grp']
        for tname in tnames:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='miRNA')
            for idx in df.index:
                if df.loc[idx, 'genes'] is None:
                    continue

                genes = set(df.loc[idx, 'genes'].split(';'))

                df.loc[idx, 'genes (org)'] = ';'.join(sorted(list(genes)))
                df.loc[idx, '#genes (org)'] = len(genes)

                genes = set.intersection(gene_names, genes)
                compl = genes - gene_names
                print(len(compl))
                df.loc[idx, 'genes'] = ';'.join(sorted(list(genes)))
                df.loc[idx, '#genes'] = len(genes)
            df.to_sql(tname, con, if_exists='replace')

    def eval(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_rna.db')
        con = sqlite3.connect(fpath)

        Y = pd.read_sql("SELECT * FROM 'Y'", con, index_col='tissue')
        X = pd.read_sql("SELECT * FROM 'X'", con, index_col='tissue')
        coeff = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='transcript_id')

        Yh = np.dot(X, coeff)
        distance = (Y - Yh).abs().values.sum() / Y.size
        print(distance)

    def compare_cline_tissue(self):
        # import matplotlib.pyplot as plt
        # from mpl_toolkits.mplot3d import Axes3D
        # fig = plt.figure()
        # ax1 = fig.add_subplot(211, projection='3d')
        # ax2 = fig.add_subplot(212, projection='3d')

        fpath_cline = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100.db')
        con_cline = sqlite3.connect(fpath_cline)
        df_cline = pd.read_sql("SELECT * FROM 'coefficient'", con_cline, index_col='transcript_id')

        fpath_rna = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_rna.db')
        con_rna = sqlite3.connect(fpath_rna)
        df_rna = pd.read_sql("SELECT * FROM 'coefficient'", con_rna, index_col='transcript_id')

        rinter = set.intersection(set(df_cline.index), set(df_rna.index))
        cinter = set.intersection(set(df_cline.columns), set(df_rna.columns))

        df_cline = df_cline.loc[rinter, cinter]
        df_rna = df_rna.loc[rinter, cinter]

        writer = pd.ExcelWriter(os.path.join(self.root, '/home/mingyu/Bioinformatics/database/Fantom/v5/cell_lines/out', 'regression_comp.xlsx'), engine='xlsxwriter')
        df_cline.to_excel(writer, sheet_name='cell_lines')
        df_rna.to_excel(writer, sheet_name='tissue')
        writer.save()
        writer.close()

        distance = ((df_cline - df_rna).abs()).values.sum() / df_rna.size
        print(distance)

        # y = np.arange(df_cline.shape[0])
        # x = np.arange(df_cline.shape[1])
        # X, Y = np.meshgrid(x, y)
        # ax1.plot_surface(X, Y, df_cline)
        # ax2.plot_surface(X, Y, df_rna)
        # plt.show()

    def distance(self):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='transcript_id')

        fpath_rna = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_rna.db')
        con_rna = sqlite3.connect(fpath_rna)
        df_rna = pd.read_sql("SELECT * FROM 'coefficient'", con_rna, index_col='transcript_id')

        mir = sorted(list(set.intersection(set(df.columns), set(df_rna.columns))))
        gene = sorted(list(set.intersection(set(df.index), set(df_rna.index))))

        df = df.loc[gene, mir]
        df_rna = df_rna.loc[gene, mir]
        diff = (df - df_rna).abs()
        distance = diff.values.sum() / df.size
        print(distance)

        xaxis = np.arange(df.size)
        plt.scatter(xaxis, df.values.flatten(), color='r', label='fantom')
        plt.scatter(xaxis, df_rna.values.flatten(), color='b', label='rna-seq')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    hostname = socket.gethostname()
    if hostname == 'mingyu-Precision-Tower-7810':
        root = '/home/mingyu/Bioinformatics'
    elif hostname == '-DLOOJR6' or hostname == '-1NLOLK4':
        root = 'D:/Bioinformatics'
    elif hostname == 'mingyu-Inspiron-7559':
        root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
    else:
        root = '/lustre/fs0/home/mcha/Bioinformatics'

    rg = Regression(root)
    if hostname == 'mingyu-Precision-Tower-7810':
        rg.to_server()
    else:
        # rg.regression(100, 'nz', 'cell_lines')
        # rg.regression_rna()
        rg.get_plugin_distance('nz')
        # for type in ['cell_lines', 'tissues']:
            # rg.regression(100, 'nz', type)
            # rg.regression_rna()
            # rg.compare()
            # rg.eval()
            # rg.compare_cline_tissue()
            # rg.report('rna')
            # rg.add_gene_name('rna')
            # rg.move()
            # rg.cross_regression(100, 'nz')

            # rg.cross_stats(100, 'nz', type=type)
            # rg.check_overlap('nz', type=type)

            # rg.cross_regression(100, 'neg')
            # rg.cross_regression(100, 'nz', N=10, type='tissues')
            # rg.cross_regression(100, 'nz', N=10, type='cell_lines')

            # rg.get_distance('nz')
            # rg.filtering(0)

        # rg.compare_tissue_cross('nz')
