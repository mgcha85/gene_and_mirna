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
        self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293']

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
            df[columns].to_sql(cline, con, if_exists='replace', index=None)

    def load_miRNA_tss(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        return pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates'", con)

    def load_gene_tss(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        df = pd.read_sql_query("SELECT * FROM 'TSS_human'", con)
        return df.groupby('chromosome')

    def main(self):
        df_mir = self.load_miRNA_tss()
        df_gene_grp = self.load_gene_tss()

        con_tag = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc.db'))
        contents = []

        N = df_mir.shape[0]
        for idx in df_mir.index:
            if idx % 100 == 0 or idx == N - 1:
                print('{:.2f}%'.format(100 * (idx + 1) / N))
            chromosome = df_mir.loc[idx, 'chromosome']
            tss = int(df_mir.loc[idx, 'tss'])
            strand = df_mir.loc[idx, 'strand']

            row = []
            for cline in self.cell_lines:
                df_tags = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND start<={tss} AND end>={tss}"
                                            "".format(cline, chromosome, tss=tss), con_tag)
                if df_tags.empty:
                    continue
                df_gene_chr = df_gene_grp.get_group(chromosome)
                midx = (df_gene_chr['start'] - tss).abs().idxmin()
                if strand == '+':
                    gene_tss = df_gene_chr.loc[midx, 'start']
                else:
                    gene_tss = df_gene_chr.loc[midx, 'end']

                df_gene_tags = pd.read_sql_query("SELECT * FROM '{}' WHERE chromosome='{}' AND NOT start>{tss} AND NOT "
                                            "end<{tss}".format(cline, chromosome, tss=gene_tss), con_tag)
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
                row = [chromosome, tss, strand, gene_tss, starts, ends, length, cline]

            if len(row) > 0:
                contents.append(row)

        df_rep = pd.DataFrame(data=contents, columns=['chromosome', 'tss', 'strand', 'gene_tss', 'tag_starts', 'tag_ends', 'length', 'cell_line'])
        df_rep.to_excel(os.path.join(self.root, 'database/Fantom/v5', self.__class__.__name__ + '.xlsx'), index=None)


if __name__ == '__main__':
    st = Search_TSSs()
    st.main()

    # fpath = os.path.join(st.root, 'database/Fantom/v5', st.__class__.__name__ + '.xlsx')
    # df = pd.read_excel(fpath)
    # df = df.dropna(subset=['chromosome'])
    # df.to_excel(fpath.replace('.xlsx', '2.xlsx'))
