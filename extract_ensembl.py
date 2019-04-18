import pandas as pd
import os
import sqlite3
import socket


class Extract_Ensembl:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def run(self):
        fpath = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export.txt')
        df = pd.read_csv(fpath, sep='\t')

        con = sqlite3.connect(fpath.replace('.txt', '.db'))
        df.to_sql('Ensembl', con, if_exists='replace', index=None)

    def run2(self):
        columns__ = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

        fpath = os.path.join(self.root, 'database/UCSC/Genes', 'genes.gtf')
        df = pd.read_csv(fpath, sep='\t', names=columns__)
        df['seqname'] = 'chr' + df['seqname'].astype(str)
        df_grp = df.groupby('source')

        con = sqlite3.connect(fpath.replace('.gtf', '.db'))
        for src, df_sub in df_grp:
            print(src)
            df_sub.to_sql(src, con, if_exists='replace', index=None)

    def compare(self):
        fpath = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export.db')
        con = sqlite3.connect(fpath)
        df_ens__ = pd.read_sql_query("SELECT * FROM 'Ensembl'", con)

        contents = []
        df_mir = pd.read_excel('miRNA.xlsx')
        for idx in df_mir.index:
            mirna = df_mir.loc[idx, 'name']
            mirna = mirna.replace('hsa', '')
            mirna = mirna.replace('-', '')

            df_ens = df_ens__[df_ens__['Gene name'] == mirna.upper()]
            if not df_ens.empty:
                mir_start = df_mir.loc[idx, 'start']
                mir_end = df_mir.loc[idx, 'end']
                chromosome = df_mir.loc[idx, 'chromosome']

                ens_tss = df_ens.iloc[0]['Transcription start site (TSS)']
                strand = df_ens.iloc[0]['Strand']

                if strand > 0:
                    mir_tss = mir_start
                else:
                    mir_tss = mir_end
                contents.append([mirna, abs(mir_tss - ens_tss)])

        df_rep = pd.DataFrame(contents, columns=['mirna', 'distance'])
        df_rep.to_excel('compare_fantom_ensembl.xlsx', index=None)


if __name__ == '__main__':
    ee = Extract_Ensembl()
    ee.compare()
