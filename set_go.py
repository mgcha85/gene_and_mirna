import pandas as pd
import sqlite3
import socket
import os
import pickle as pkl
from time import sleep
from bs4 import BeautifulSoup
from Database import Database


class set_go:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def set_input(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out')
        fpath = os.path.join(dirname, 'regression_100.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'result'", con, index_col='miRNA')
        for i, mir in enumerate(df.index):
            print(i + 1, df.shape[0])
            genes = df.loc[mir, 'gene_name']
            genes = genes.split(';')
            with open(os.path.join(dirname, 'gene_set', mir + '.txt'), 'wt') as f:
                f.write('\n'.join(genes))

    def load_gene_list(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/gene_set')
        flist = os.listdir(dirname)

        def get_flist(flist, dirname):
            for fname in flist:
                with open(os.path.join(dirname, fname), 'rt') as f:
                    yield os.path.splitext(fname)[0], f.read()
        return get_flist(flist, dirname)

    def get_bg_genes(self):
        fpath = os.path.join(self.root, 'database/gencode', 'high_correlated_genes.txt')
        with open(fpath, 'rt') as f:
            return f.read()

    def submit_data(self):
        from mechanize import Browser
        br = Browser()

        bg_genes = self.get_bg_genes()
        gene_list = self.load_gene_list()

        responses = {}
        for i, glist in enumerate(gene_list):
            print(i+1, glist[0])
            br.open("http://cbl-gorilla.cs.technion.ac.il/")
            for f in br.forms():
                if f.attrs['name'] == 'gorilla':
                    br.form = f

            target_set = br.form.find_control(name="target_set")
            target_set.value = glist[1]
            background_set = br.form.find_control(name="background_set")
            background_set.value = bg_genes

            run_mode = br.form.find_control(name="run_mode")
            run_mode.items[0].selected = False
            run_mode.items[1].selected = True

            ntry, n = 0, 3
            while ntry < n:
                try:
                    response = br.submit(type='submit')
                    break
                except:
                    print('connection error')
                    print('[{} / {}] try'.format(ntry+1, n))
                    sleep(1)
                    ntry += 1
                    continue

            if ntry == n:
                print('{} failed'.format(glist[0]))
                continue

            res = response.read().decode('utf-8')
            with open(os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', glist[0] + '.html'), 'wt') as f:
                f.write(res)
            sleep(0.5)

    def extract_genes(self):
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

        dirname = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result')
        flist = os.listdir(dirname)
        con = sqlite3.connect(os.path.join(dirname, 'gene_list.db'))
        for fname in flist:
            if fname.endswith('.html'):
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
                    df_rep.to_sql(mir, con, index=None, if_exists='replace')
                else:
                    print('{} does not have target genes'.format(mir))

    def summary(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', 'gene_list.db')
        con = sqlite3.connect(fpath)
        tlist = Database.load_tableList(con)

        contents = []
        for tname in tlist:
            df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
            contents.append([tname, df.iloc[0]['Genes'], df.iloc[0]['P-value']])

        df_res = pd.DataFrame(contents, columns=['miRNA', 'genes', 'p-value'])

        fpath = os.path.join(self.root, 'database', 'target_genes.db')
        df_res.to_sql('go_result', sqlite3.connect(fpath), index=None)

    def gene_to_tr(self):
        fpath = os.path.join(self.root, 'database/gencode', 'gencode.v32lift37.annotation_attr.db')
        con = sqlite3.connect(fpath)

        tlist = Database.load_tableList(con)
        dfs = []
        for tname in tlist:
            dfs.append(pd.read_sql("SELECT transcript_id, gene_name FROM '{}' WHERE transcript_id IS NOT NULL "
                                   "AND feature='transcript'".format(tname), con, index_col='gene_name'))
        return pd.concat(dfs)

    def add_corr(self):
        def gene_to_tr(df, genes, con):
            trs = []
            for gene in genes:
                if gene not in df.index:
                    continue
            
                df_tr = df.loc[gene]
                sql = ["transcript_id='{}'".format(x) for x in df_tr['transcript_id']]
                sql = ' OR '.join(sql)
                sql = "SELECT * FROM 'corr' WHERE " + sql
                corr_coef = pd.read_sql(sql, con, index_col='transcript_id')
                if corr_coef.empty:
                    continue

                mean = corr_coef.mean(axis=1)
                trs.append(mean.idxmax())
            return trs

        fpath_ref = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100.db')
        con_ref = sqlite3.connect(fpath_ref)

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out/go_result', 'gene_list_out.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql("SELECT * FROM 'go_result'", con, index_col='miRNA')
        genes = df['genes'].str.split(';')

        df_tr = self.gene_to_tr()
        for i, mir in enumerate(df.index):
            print('{:,d} / {:,d}'.format(i+1, df.shape[0]))
            gname = genes[mir]
            tname = gene_to_tr(df_tr, gname, con_ref)

            sql = ["transcript_id='{}'".format(x) for x in tname]
            sql = ' OR '.join(sql)
            sql = "SELECT * FROM 'corr' WHERE " + sql
            corr_coef = pd.read_sql(sql, con_ref, index_col='transcript_id')
            if corr_coef.empty:
                continue
            df.loc[mir, 'corr'] = ';'.join(corr_coef[mir].round(4).astype(str).values)
            df.loc[mir, 'corr_stats'] = '{:0.4f};{:0.4f}'.format(corr_coef[mir].abs().mean(), corr_coef[mir].abs().std())
        df.to_sql('go_result', con, if_exists='replace')


if __name__ == '__main__':
    sg = set_go()
    # sg.submit_data()
    # sg.extract_genes()
    sg.summary()
