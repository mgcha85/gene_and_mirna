import pandas as pd
import socket
import os
import sqlite3



class Search_TSSs:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        elif 'evc' in hostname:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        else:
            print('wrong option')
            return
        self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293', 'MCF7']

    def split_table_by_cell_lines(self):
        fname__ = 'hg19.cage_peak_phase1and2combined_counts.osc'
        dirname = os.path.join(self.root, 'database/Fantom/v5')
        fpath = os.path.join(dirname, fname__ + '.csv')
        df = pd.read_csv(fpath)
        # columns = list(df.columns)
        # with open('columns.txt', 'wt') as f:
        #     f.write('\n'.join(columns))
        # exit(1)

        con = sqlite3.connect(os.path.join(dirname, 'hg19_cage_peak_phase1and2combined_counts_osc.db'))
        for cline in self.cell_lines:
            print(cline)
            columns = ['chromosome', 'start', 'end', 'strand']
            for col in df.columns:
                if cline in col:
                    columns.append(col)
            df = df.drop_duplicates()
            df[columns].to_sql(cline, con, if_exists='replace', index=None)

    def load_miRNA_tss(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        return pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE Type='intronic'", con)

    def load_gene_tss(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        return pd.read_sql_query("SELECT * FROM 'human_gene'", con, index_col='gene')

    def main(self):
        df_mir = self.load_miRNA_tss()
        df_gene_grp = self.load_gene_tss()

        con_tag = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc.db'))
        contents = []
        columns = ['chromosome', 'pre_start', 'pre_end', 'tss', 'strand', 'gene_tss', 'gene_start', 'gene_end',
                   'tag_starts', 'tag_ends', 'length', 'gene', 'cell_line']

        N = df_mir.shape[0]
        for idx in df_mir.index:
            if idx % 100 == 0 or idx == N - 1:
                print('{:.2f}%'.format(100 * (idx + 1) / N))
            chromosome = df_mir.loc[idx, 'chromosome']
            tss = int(df_mir.loc[idx, 'tss'])
            strand = df_mir.loc[idx, 'strand']
            genes = df_mir.loc[idx, 'Gene']
            pre_start = df_mir.loc[idx, 'start']
            pre_end = df_mir.loc[idx, 'end']
            genes = genes.split(',')

            row = []
            for gene in genes:
                for cline in self.cell_lines:
                    df_tags = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND start<={tss} AND end>={tss}"
                                                "".format(cline, chromosome, tss=tss), con_tag)
                    if df_tags.empty:
                        continue
                    try:
                        df_gene_name = df_gene_grp.loc[gene]
                    except Exception as e:
                        print(e)
                        continue

                    gene_start = df_gene_name['start']
                    gene_end = df_gene_name['end']
                    if strand == '+':
                        gene_tss = df_gene_name['start']
                    else:
                        gene_tss = df_gene_name['end']

                    df_gene_tags = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND NOT start>{tss} AND NOT "
                                                "end<{tss}".format(cline, chromosome, tss=gene_tss), con_tag)
                    if df_gene_tags.empty:
                        continue
                    df_gene_tags = df_gene_tags[df_gene_tags.iloc[:, 4:].sum(axis=1) > 0]
                    if df_gene_tags.empty:
                        continue

                    tsss = sorted([gene_tss, tss])
                    start = tsss[0]
                    end = tsss[1]
                    length = end - start
                    df_tags = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND NOT start<={} AND NOT "
                                                "end>={}".format(cline, chromosome, start, end), con_tag)

                    starts = ';'.join(list(df_tags['start'].values.astype(str)))
                    ends = ';'.join(list(df_tags['end'].values.astype(str)))
                    row.append([chromosome, pre_start, pre_end, tss, strand, gene_tss, gene_start, gene_end, starts,
                                ends, length, gene, cline])

            if len(row) > 0:
                contents.extend(row)

        df_rep = pd.DataFrame(data=contents, columns=columns)
        df_rep.to_excel(os.path.join(self.root, 'database/Fantom/v5', self.__class__.__name__ + '.xlsx'), index=None)


if __name__ == '__main__':
    st = Search_TSSs()
    st.main()

    # fpath = os.path.join(st.root, 'database/Fantom/v5', st.__class__.__name__ + '.xlsx')
    # df = pd.read_excel(fpath)
    # df = df.dropna(subset=['chromosome'])
    # df.to_excel(fpath.replace('.xlsx', '2.xlsx'))
