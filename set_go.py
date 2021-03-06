import pandas as pd
import sqlite3
import socket
import os
import pickle as pkl
from time import sleep
from bs4 import BeautifulSoup
from Database import Database
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom


class set_go:
    def __init__(self, root):
        self.root = root

    def to_server(self, root, tag):
        from Server import Server
        import sys

        which = 'stokes'
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

    def set_input(self, hbw, opt):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out')
        fpath = os.path.join(dirname, 'regression_{}_{}.db'.format(hbw, opt))

        # dirname = os.path.join(self.root, 'database')
        # fpath = os.path.join(dirname, 'target_genes.db')

        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'result'", con, index_col='miRNA')
        # df = pd.read_sql("SELECT * FROM 'target_scan_grp'", con, index_col='miRNA')
        # df = pd.read_sql("SELECT * FROM 'Amlan_result'", con, index_col='miRNA')

        con_out = sqlite3.connect(os.path.join(dirname, 'gene_set', 'db_{}_{}.sqlite'.format(hbw, opt)))
        # con_out = sqlite3.connect(os.path.join(dirname, 'TargetScan/go_result', 'db_org.sqlite'))

        df['gene_name'] = df['gene_name'].str.split(';')
        df['Transcripts'] = df['Transcripts'].str.split(';')

        hgenes = self.get_bg_genes().split('\n')
        # hgenes = self.get_bg_genes()
        # df['gene_name'] = df['gene_name'].str.split(';')
        # df['Transcripts'] = df['Transcripts'].str.split(';')

        for i, mir in enumerate(df.index):
            print(i + 1, df.shape[0])
            genes = df.loc[mir, 'gene_name']
            trans = df.loc[mir, 'Transcripts']

            df_res = pd.DataFrame(columns=['genes', 'transcripts'])
            df_res['genes'] = genes
            df_res['transcripts'] = trans

            df_res = df_res.set_index('genes', drop=True)
            genes = set.intersection(set(df_res.index), set(hgenes))
            df_res = df_res.loc[genes]
            df_res.to_sql(mir, con_out, if_exists='replace')

    def get_bg_genes(self):
        fpath = os.path.join(self.root, 'database/gencode', 'high_correlated_fan_rna_100.xlsx')
        df = pd.read_excel(fpath)
        return '\n'.join(sorted(list(set(df['gene_name']))))
        # with open(fpath, 'rt') as f:
        #     return f.read()

    def get_con_mirna(self):
        fpath = os.path.join(self.root, 'database', 'consistent_miRNA_330.db')
        con = sqlite3.connect(fpath)
        return pd.read_sql("SELECT DISTINCT miRNA FROM 'mean_tss'", con)

    def get_bg_genes_each(self):
        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'target_scan'", con)
        return '\n'.join(df.drop_duplicates(subset=['Gene Symbol'])['Gene Symbol'])

    def submit_data(self, hbw, opt, bg=True, mode='proc'):
        from mechanize import Browser
        br = Browser()

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/gene_set', 'db_{}_{}.sqlite'.format(hbw, opt))
        # fpath = os.path.join(self.root, 'database/mirTarBase/go_result', 'db_org.sqlite')
        con = sqlite3.connect(fpath)
        mlist = Database.load_tableList(con)

        con_mirs = self.get_con_mirna()
        mlist = sorted(list(set.intersection(set(con_mirs['miRNA']), set(mlist))))

        fails = []
        for i, mir in enumerate(mlist):
            print(i+1, len(mlist))

            genes = pd.read_sql("SELECT DISTINCT genes FROM '{}'".format(mir), con)
            if genes.empty:
                continue

            br.open("http://cbl-gorilla.cs.technion.ac.il/")
            for f in br.forms():
                if f.attrs['name'] == 'gorilla':
                    br.form = f

            target_set = br.form.find_control(name="target_set")
            target_set.value = '\n'.join(genes['genes'])
            run_mode = br.form.find_control(name="run_mode")
            if bg:
                run_mode.items[0].selected = False
                run_mode.items[1].selected = True
                background_set = br.form.find_control(name="background_set")
                bg_genes = self.get_bg_genes()
                background_set.value = bg_genes

            db = br.form.find_control(name="db")
            for i in range(4):
                if i == mode:
                    db.items[i].selected = True
                else:
                    db.items[i].selected = False

            ntry, n = 0, 3
            while ntry < n:
                try:
                    response = br.submit(type='submit')
                    break
                except:
                    print('connection error')
                    print('[{} / {}] try'.format(ntry+1, n))
                    sleep(0.5)
                    ntry += 1
                    continue

            if ntry == n:
                print('{} failed'.format(mir))
                fails.append(mir)
                continue

            res = response.read().decode('utf-8')
            dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', opt, db.items[mode].name)
            if not os.path.exists(dirname):
                os.mkdir(dirname)

            with open(os.path.join(dirname, mir + '.html'), 'wt') as f:
                f.write(res)
            sleep(1)

        if fails:
            with open('failure.txt', 'wt') as f:
                f.write('\n'.join(fails))

    def extract_genes(self, hbw, opt, mode):
        def trim_genes(df):
            for idx in df.index:
                genes = df.loc[idx, 'Genes']
                genes = genes.split('Show genes')[1]
                genes = genes.split('\n')

                gl = []
                for gene in genes:
                    if len(gene) < 2:
                        continue
                    gl.append(gene.split(' - ')[0].strip())
                df.loc[idx, 'Genes'] = ';'.join(gl)
            return df

        # thres = 0.1
        # dirname = os.path.join(self.root, 'database/mirTarBase/go_result')
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', opt, mode)
        flist = [x for x in os.listdir(dirname) if x.endswith('.html')]
        con = sqlite3.connect(os.path.join(dirname, 'gene_list_{}_{}.db'.format(hbw, opt)))
        for fname in flist:
            fpath = os.path.join(dirname, fname)

            with open(fpath, 'rt') as f:
                page = f.read()
            soup = BeautifulSoup(page, "lxml")
            script = soup.findAll("tr")

            table = []
            for row in script[1:]:
                ele = row.findAll("td")
                row_ele = []
                for e in ele:
                    row_ele.append(e.text.strip())
                table.append(row_ele)

            mir = os.path.splitext(fname)[0]
            if len(table) > 0:
                df_rep = pd.DataFrame(table)
                df_rep.columns = df_rep.iloc[0]
                df_rep = df_rep[1:]
                df_rep = trim_genes(df_rep)
                df_rep[['P-value', 'FDR q-value']] = df_rep[['P-value', 'FDR q-value']].astype(float)
                # df_rep = df_rep[df_rep['FDR q-value'] < thres]
                df_rep = df_rep.rename(columns={'FDR q-value': 'q-value'})
                if df_rep.empty:
                    print('{} does not have target genes'.format(mir))
                    continue
                df_rep.sort_values('q-value').to_sql(mir, con, index=None, if_exists='replace')
            else:
                print('{} does not have target genes'.format(mir))

    def summary(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', 'gene_list.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        contents = []
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}' LIMIT 1".format(tname), con)
            contents.append([tname, df.iloc[0]['Genes'], df.iloc[0]['q-value']])

        df_res = pd.DataFrame(contents, columns=['miRNA', 'genes', 'q-value'])

        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        df_res.to_sql('go_result', sqlite3.connect(fpath), index=None, if_exists='replace')

    def gene_to_tr(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v32lift37.annotation_attr.db')
        con = sqlite3.connect(fpath)

        tlist = Database.load_tableList(con)
        dfs = []
        for tname in tlist:
            dfs.append(pd.read_sql("SELECT transcript_id, gene_name FROM '{}' WHERE transcript_id IS NOT NULL "
                                   "AND feature='transcript'".format(tname), con, index_col='gene_name'))
        return pd.concat(dfs)

    def add_corr(self, hbw=100):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines')
        con_corr = sqlite3.connect(os.path.join(dirname, 'out', 'regression_{}_{}.db'.format(hbw, opt)))

        fpath_tid = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/gene_set', 'db.sqlite')
        con_tid = sqlite3.connect(fpath_tid)

        fpath = os.path.join(dirname, 'out/go_result', 'gene_list.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)
        for mir in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(mir), con)
            df_tid = pd.read_sql("SELECT * FROM '{}'".format(mir), con_tid, index_col='genes')

            for i, genes in enumerate(df['Genes'].str.split(';')):
                corr = np.zeros(len(genes))
                for j, gene in enumerate(genes):
                    tid = df_tid.loc[gene, 'transcripts']
                    if not isinstance(tid, str):
                        tid = tid.iloc[0]
                    df_corr = pd.read_sql("SELECT * FROM 'corr' WHERE transcript_id='{tid}'".format(mir=mir, tid=tid), con_corr, index_col='transcript_id')
                    corr[j] = df_corr.loc[tid, mir]

                corr = corr.round(4)
                df.loc[i, 'corr'] = ';'.join(corr.astype(str))
                stats = np.array([corr.mean(), corr.std(), np.median(corr), corr.min(), corr.max()])
                df.loc[i, 'stats'] = ';'.join(stats.round(4).astype(str))

            df.to_sql(mir, con, if_exists='replace', index=None)

    def confusion(self):
        report = []
        for dirname in ['mirTarBase', 'TargetScan', 'Amlan', 'Lasso']:
            fpath = os.path.join(self.root, 'database', dirname, 'go_result', 'gene_list.db')
            con = sqlite3.connect(fpath)

            contents = []
            for tname in Database.load_tableList(con):
                df = pd.read_sql("SELECT * FROM '{}' LIMIT 1".format(tname), con)
                contents.append([tname, df.loc[0, 'q-value']])
            report.append(pd.DataFrame(contents, columns=['miRNA', 'q-value ({})'.format(dirname)]).set_index('miRNA'))
        df_res = pd.concat(report, axis=1)
        df_res.to_excel(os.path.join(self.root, 'database', 'confusion.xlsx'))

    def result(self, hbw, opt, mode):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}_{}.db'.format(hbw, opt))
        con = sqlite3.connect(fpath)

        con_go = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', opt, mode, 'gene_list_{}_{}.db'.format(hbw, opt)))
        df = pd.read_sql("SELECT * FROM 'result'", con, index_col='miRNA')
        result = []
        for mir in df.index:
            if Database.checkTableExists(con_go, mir):
                df_go = pd.read_sql("SELECT * FROM '{}' LIMIT 1".format(mir), con_go)
                go_gene = df_go.loc[0, 'Genes']
            else:
                go_gene = None
            result.append([mir, df.loc[mir, 'Transcripts'], df.loc[mir, 'gene_name'], go_gene, df_go.loc[0, 'q-value']])
        pd.DataFrame(result, columns=['miRNA', 'Transcripts', 'gene_lasso', 'gene_go', 'q-value']).to_sql('lasso_go_{}'.format(mode), con, index=False, if_exists='replace')

    def merge_results(self, hbw, opt):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}_{}.db'.format(hbw, opt))
        con = sqlite3.connect(fpath)

        dfs = []
        sig_mir = []
        for tname in ['lasso_go_comp', 'lasso_go_func', 'lasso_go_proc']:
            cname = 'q-value ({})'.format(tname.split('_')[-1])
            df = pd.read_sql('SELECT "miRNA", "q-value" AS "{}" FROM "{}"'.format(cname, tname),
                        con, index_col='miRNA')
            sig_mir.append(set(df[df[cname] < 0.01].index))
            dfs.append(df)
        df_res = pd.concat(dfs, axis=1)
        df_res['significant'] = 'N'
        sig_mir = set.union(*sig_mir)
        df_res.loc[sig_mir, 'significant'] = 'Y'
        df_res.to_sql('lasso_go', con, if_exists='replace')

    def add_mir_type(self):
        fpath = os.path.join(self.root, 'database', 'consistent_miRNA_330.db')
        con = sqlite3.connect(fpath)
        con_pred = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'predictions_processed_result.db'))

        df = pd.read_sql("SELECT miRNA FROM 'mean_tss'", con)
        for tname in Database.load_tableList(con_pred):
            df_pred = pd.read_sql("SELECT * FROM '{}'".format(tname), con_pred, index_col='miRNA')
            con_out = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'predictions_processed_result2.db'))
            mirna = set.intersection(set(df['miRNA']), set(df_pred.index))
            df_pred.loc[mirna, 'feature'] = 'precursor'
            idx = df_pred[df_pred['feature'] != 'precursor'].index
            df_pred.loc[idx, 'feature'] = 'miRNA'
            df_pred.to_sql(tname, con_out, if_exists='replace')

    def to_pre_mir(self):
        fpath = os.path.join(self.root, 'database/mirbase/22', 'mirbase.db')
        con = sqlite3.connect(fpath)
        con_pred = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'predictions_processed_result.db'))
        con_out = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'predictions_processed_out.db'))

        df = pd.read_sql("SELECT * FROM 'hsa_liftover_hg19'", con, index_col='ID')
        for tname in Database.load_tableList(con_pred):
            if '_pre' in tname:
                df_pred = pd.read_sql("SELECT * FROM '{}'".format(tname), con_pred)
                for idx in df_pred.index:
                    mir = df_pred.loc[idx, 'pre-miRNA']

                    id = df[df['Name'] == mir].index
                    if len(id) == 0:
                        continue
                    else:
                        id = id[0]

                    derives = df.loc[id, 'Derives_from']
                    if derives:
                        name = df.loc[derives, 'Name']
                        df_pred.loc[idx, 'miRNA'] = name
                df_pred.to_sql(tname, con_out, index=None, if_exists='replace')

    def move(self, hbw, opt):
        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed.db')
        con = sqlite3.connect(fpath)
        fpath_out = os.path.join(self.root, 'database/target_genes', 'predictions_processed_{}_{}.db'.format(hbw, opt))
        con_out = sqlite3.connect(fpath_out)
        for tname in Database.load_tableList(con):
            if '_pre' in tname and tname != 'go_pre':
                df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
                df.to_sql(tname, con_out, index=False, if_exists='replace')

    def hyper_test(self, hbw, opt):
        high_genes = self.get_bg_genes().split('\n')
        out_path = os.path.join(self.root, 'database/target_genes', 'hyper_test_lasso_{}_{}.xlsx'.format(hbw, opt))

        con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}_{}.db'.format(hbw, opt)))
        con_pred = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'predictions_processed_{}_{}.db'.format(hbw, opt)))
        N = len(high_genes)

        df = pd.read_sql("SELECT miRNA, gene_lasso FROM 'lasso_go'", con, index_col='miRNA')
        thres = 1 / df.shape[0]
        tlist = [x for x in Database.load_tableList(con_pred) if '_pre' in x]
        for tname in tlist:
            print(tname)
            df_pred = pd.read_sql("SELECT * FROM '{}'".format(tname), con_pred, index_col='pre-miRNA')
            mir = sorted(list(set.intersection(set(df.index), set(df_pred.index))))
            df.loc[mir, 'genes ({})'.format(tname)] = df_pred.loc[mir, 'genes']

        cols = [col for col in df.columns if 'genes (' in col]
        for mir in df.index:
            gene_go = df.loc[mir, 'gene_lasso']
            if gene_go is None:
                continue
            gene_lasso = set(gene_go.split(';'))

            target_genes = {}
            for col in cols:
                tg = df.loc[mir, col]
                key = col.replace('genes (', '')[:-1]
                if isinstance(tg, float):
                    target_genes[key] = set([])
                else:
                    target_genes[key] = set.intersection(set(high_genes), set(tg.split(';')))

            for vname, method in target_genes.items():
                vname = vname.replace('_pre', '')

                if vname == 'common':
                    df_add = pd.read_sql("SELECT * FROM '{}'".format(vname + '_pre'), con_pred, index_col='pre-miRNA')
                    idx = set.intersection(set(df.index), set(df_add.index))
                    df.loc[idx, 'source (common)'] = df_add.loc[idx, 'source']

                q = len(set.intersection(gene_lasso, method))
                m = len(method)
                n = N - m
                k = len(gene_lasso)
                p = hypergeom.sf(q-1, n+m, m, k)
                df.loc[mir, 'p-val ({})'.format(vname)] = None
                if q > 0:
                    df.loc[mir, 'p-val ({})'.format(vname)] = 1 - p
                df.loc[mir, 'ele ({})'.format(vname)] = ','.join(map(str, [q, m, n, k]))

        df = self.add_go(df, hbw, opt)
        df['# significant'] = 0
        values = df[[col for col in df.columns if 'p-val (' in col]].values
        idx = np.where(values < thres)
        if len(idx) == 0:
            df.drop(cols + ['gene_lasso'], axis=1).to_excel(out_path)
        else:
            add_col = {}
            for i in idx[0]:
                add_col[i] = add_col.get(i, 0) + 1
            for r, v in add_col.items():
                df.iloc[r, -1] = v

            df.sort_index().drop(cols + ['gene_lasso'], axis=1).to_excel(out_path)

    def hyper_gsea(self):
        def get_gene_set(dirname):
            flist = [x for x in os.listdir(dirname) if x.endswith('.xls') and 'HALLMARK' in x]

            contents = []
            for fname in flist:
                fpath = os.path.join(dirname, fname)
                df = pd.read_csv(fpath, sep='\t')
                contents.append(';'.join(df['PROBE']))
            return contents

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'lasso_go'", con, index_col='miRNA')

        # tr_id to gene_name
        fpath = os.path.join(self.root, 'database/gencode', 'gene_tid.csv')
        df_list = pd.read_csv(fpath, index_col=0, sep='\t')
        bg_gene = pd.read_sql("SELECT * FROM 'X'", con, index_col='tid')
        bg_gene = df_list.loc[bg_gene.index, 'gene_name'].dropna()
        bg_gene = set(bg_gene)

        dirname = '/home/mingyu/gsea_home/output/feb13/my_analysis.Gsea.1581949269965'
        gene_sets = get_gene_set(dirname)
        for mir in df.index:
            genes = df.loc[mir, 'gene_lasso'].split(';')
            for i, gs in enumerate(gene_sets):
                gene_set = set(gs.split(';'))
                target_genes = set(genes)
                gene_inter = set.intersection(gene_set, target_genes)

                q = len(gene_inter)
                if q <= 0:
                    continue

                m = len(gene_set)
                n = len(bg_gene)
                k = len(target_genes)
                p = hypergeom.sf(q-1, n+m, m, k)
                df.loc[mir, 'p-val (GSEA set{})'.format(i)] = 1 - p

        data = df.iloc[:, 4:].values
        thres = 0.01 / df.shape[0]
        rix, cix = np.where(data < thres)

        sets = {}
        for r, c in zip(rix, cix):
            if r in sets:
                sets[r].append(str(c+1))
            else:
                sets[r] = [str(c+1)]

        contents = []
        for r, c in sets.items():
            contents.append([df.index[r], ';'.join(c)])
        df_add = pd.DataFrame(contents, columns=['tid', 'set_num']).set_index('tid')

        df['significant'] = 0
        df['set_num'] = None
        df.loc[df_add.index, 'set_num'] = df_add.loc[:, 'set_num']
        df.loc[df_add.index, 'significant'] = 1
        df.to_sql('lasso_gsea', con, if_exists='replace')

    def plot_hyper_test(self, hbw, opt):
        fpath = os.path.join(self.root, 'database/target_genes', 'hyper_test_lasso_{}_{}.xlsx'.format(hbw, opt))
        df = pd.read_excel(fpath, index_col=0)
        columns = [x for x in df.columns if 'p-val' in x]

        xaxis = np.arange(df.shape[0])
        for i, col in enumerate(columns):
            yaxis = df[col]
            label = col.replace('p-val', '').replace('(', '').replace(')', '')
            plt.subplot(4, 2, i+1)
            plt.scatter(xaxis, yaxis)
            plt.title(label)
        plt.xticks(xaxis[::10], df.index[::10], rotation=30, fontsize=5)
        # plt.legend()
        plt.show()

    def add_go(self, df__, hbw, opt):
        con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', opt, 'gene_list_{}_{}.db'.format(hbw, opt)))

        contents = []
        for tname in Database.load_tableList(con):
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            contents.append([tname, df['q-value'].iloc[0]])
        df = pd.DataFrame(data=contents, columns=['pre-miRNA', 'q-value']).set_index('pre-miRNA')

        idx = set.intersection(set(df.index), set(df__.index))
        df__.loc[idx, 'q-value (GO)'] = df.loc[idx, 'q-value']
        df__ = df__.sort_values(by='q-value (GO)')
        df__['q significant'] = df__['q-value (GO)'] * np.arange(1, df__.shape[0]+1)
        return df__

    def to_tg(self, hbw, opt, mode):
        con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_{}_{}.db'.format(hbw, opt)))
        df = pd.read_sql("SELECT miRNA, gene_go FROM 'lasso_go' WHERE gene_go IS NOT NULL", con, index_col='miRNA')
        df.index.name = "pre-miRNA"
        df.columns = ["genes"]

        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed_result_{}_{}.db'.format(hbw, opt))
        out_con = sqlite3.connect(fpath)
        df.to_sql('go_pre', out_con, if_exists='replace')

    def figure(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'hyper_test_lasso.xlsx')
        df = pd.read_excel(fpath, index_col=0)
        xaxis = range(df.shape[0])

        plt.figure(figsize=(12, 8))
        ttl = []
        for vname in ['miranda', 'rna22', 'ts', 'common']:
            plt.scatter(xaxis, df['p-val ({})'.format(vname)], label=vname)
            ttl.append('[{}] mean: {:0.2f}, std: {:0.2f}'.format(vname, df['p-val ({})'.format(vname)].mean(), df['p-val ({})'.format(vname)].std()))

        plt.xticks(xaxis[::10], df.index[::10], rotation=30, fontsize=6)
        plt.legend()
        plt.grid()
        plt.title('\n'.join(ttl))
        plt.savefig(os.path.join(self.root, 'database/target_genes', 'hyper_test_lasso.png'))

    def add_pre_mir(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed.db')
        con = sqlite3.connect(fpath)

        con_mir = sqlite3.connect(os.path.join(self.root, 'database/target_genes', 'mirna_primirna_mapping.db'))

        # tlist = Database.load_tableList(con)
        tlist = ['miranda_grp_mir', 'rna22_grp_mir', 'ts_grp_mir', 'union']
        for tname in tlist:
            print(tname)
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='mir')
            for mir in df.index:
                mid = pd.read_sql("SELECT * FROM 'mirna_primirna_mapping' WHERE mirna_name='{}'".format(mir), con_mir)
                if mid.empty:
                    print(mir)
                    continue
                df.loc[mir, 'pre-miRNA'] = ';'.join(mid['primirna_name'])
            df.to_sql(tname, con, if_exists='replace')

    def convert_pre(self):
        fpath = os.path.join(self.root, 'database/target_genes', 'predictions_processed.db')
        con = sqlite3.connect(fpath)
        # tlist = Database.load_tableList(con)

        tlist = ['intersection', 'miranda_grp_mir', 'rna22_grp_mir', 'ts_grp_mir', 'union']
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            contents = []
            for pmir, df_grp in df.groupby('pre-miRNA'):
                if ';' in pmir:
                    continue

                genes = ';'.join(df_grp['genes'])
                genes = genes.split(';')
                genes = set(genes)
                contents.append([pmir, ';'.join(genes)])
            df_res = pd.DataFrame(contents, columns=['pre-miRNA', 'genes'])
            df_res.to_sql(tname + '_pre', con, index=False, if_exists='replace')

    def validation(self):
        # import matplotlib.pyplot as plt
        # from mpl_toolkits.mplot3d import Axes3D
        # fig = plt.figure()
        # ax1 = fig.add_subplot(211, projection='3d')
        # ax2 = fig.add_subplot(212, projection='3d')

        # 1
        # prediction RNA-seq Y using FANTOM coefficient
        fpath_fan = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz.db')
        con_fan = sqlite3.connect(fpath_fan)
        df_fan1 = pd.read_sql("SELECT * FROM 'coefficient'", con_fan, index_col='transcript_id')
        # df_fanx1 = pd.read_sql("SELECT * FROM 'X'", con_fan1, index_col='transcript_id')
        # df_fany1 = pd.read_sql("SELECT * FROM 'Y'", con_fan1, index_col='miRNA')

        # 2
        # prediction RNA-seq Y using FANTOM coefficient
        df_fan = pd.read_sql("SELECT * FROM 'coefficient'", con_fan, index_col='transcript_id')
        # df_fanx = pd.read_sql("SELECT * FROM 'X'", con_fan, index_col='transcript_id')
        # df_fany = pd.read_sql("SELECT * FROM 'Y'", con_fan, index_col='miRNA')

        fpath_rna = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_rna.db')
        con_rna = sqlite3.connect(fpath_rna)
        df_rna = pd.read_sql("SELECT * FROM 'coefficient'", con_rna, index_col='transcript_id')
        trs = set.intersection(set(df_fan.index), set(df_rna.index))

        df_fan1 = df_fan1.loc[trs, df_rna.columns]
        df_rnax = pd.read_sql("SELECT * FROM 'X'", con_rna, index_col='transcript_id')
        df_rnax = df_rnax.loc[trs]
        df_rnay = pd.read_sql("SELECT * FROM 'Y'", con_rna, index_col='Name')

        yh = df_fan1.T.dot(df_rnax)
        distance = (df_rnay - yh).abs().values.sum() / df_rnay.size
        print('1. prediction RNA-seq Y using FANTOM coefficient: {:0.6f}'.format(distance))

        df_fan = df_fan.loc[trs, df_rna.columns]
        yh = df_fan.T.dot(df_rnax)
        distance = (df_rnay - yh).abs().values.sum() / df_rnay.size
        print('2. prediction RNA-seq Y using FANTOM coefficient: {:0.6f}'.format(distance))

        # y, x = df_rnay.shape
        # X, Y = np.meshgrid(range(x), range(y))
        # ax1.plot_surface(X, Y, df_rnay)
        # ax2.plot_surface(X, Y, yh)
        # plt.show()

        # distance between RNA-seq actual Y and predicted Y
        df_rna = df_rna.loc[trs]
        yh = df_rna.T.dot(df_rnax)
        distance = (df_rnay - yh).abs().values.sum() / df_rnay.size
        print('3. distance between RNA-seq actual Y and predicted Y: {:0.6f}'.format(distance))


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

    sg = set_go(root)
    if hostname == 'DESKTOP-DLOOJR6' or hostname == '-1NLOLK4':
        sg.to_server(root, "")
    else:
        # sg.set_input(100, 'neg')
        # sg.submit_data(100, 'neg', bg=True)
        # sg.extract_genes(100, 'neg')

        # sg.result(100, 'neg')
        # sg.to_tg(100, 'neg')

        # sg.hyper_gsea()
        sg.add_pre_mir()
        sg.convert_pre()
        sg.move(100, 'nz')
        sg.hyper_test(100, 'nz')
        # sg.plot_hyper_test(100, 'neg')

        # sg.figure()
        # sg.validation()
        # sg.confusion()
        # sg.add_corr()
        # sg.summary()
