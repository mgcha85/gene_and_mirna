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
        self.cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293', 'MCF-7']

    def split_table_by_cell_lines(self, cell_lines):
        fname__ = 'hg19.cage_peak_phase1and2combined_counts.osc'
        dirname = os.path.join(self.root, 'database/Fantom/v5')
        fpath = os.path.join(dirname, fname__ + '.csv')
        df = pd.read_csv(fpath)
        # columns = list(df.columns)
        # with open('columns.txt', 'wt') as f:
        #     f.write('\n'.join(columns))
        # exit(1)

        con = sqlite3.connect(os.path.join(dirname, 'hg19_cage_peak_phase1and2combined_counts_osc.db'))
        for cline in cell_lines:
            print(cline)
            columns = ['chromosome', 'start', 'end', 'strand']
            for col in df.columns:
                if cline in col:
                    columns.append(col)

            df[columns].to_sql(cline, con)

    def load_miRNA_tss(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        return pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates'", con)

    def load_gene_tss(self):
        fpath = os.path.join(self.root, 'database', 'fantom5.db')
        con = sqlite3.connect(fpath)
        return pd.read_sql_query("SELECT * FROM 'TSS_human'", con)

    def main(self):
        dirname = os.path.join(self.root, 'database/Fantom/v5')
        fname__ = 'hg19.cage_peak_phase1and2combined_counts.osc'
        df_mir = self.load_gene_tss()
        df_gene = self.load_gene_tss()

        for idx in df_mir.index:
            chromosome = df_mir.loc[idx, 'chromosome']
            tss = df_mir.loc[idx, 'tss']
            strand = df_mir.loc[idx, 'strand']

        for cline in self.cell_lines:
            fpath = os.path.join(dirname, fname__.format(cline))
            df = pd.read_csv(fpath)




if __name__ == '__main__':
    st = Search_TSSs()
    st.main()
