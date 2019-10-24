import pandas as pd
import os
import sqlite3
import socket


class Extract_Ensembl:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
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
        df_mir = pd.read_excel('gene.xlsx')
        for idx in df_mir.index:
            print('{} / {}'.format(idx + 1, df_mir.shape[0]))
            gene = df_mir.loc[idx, 'name']

            df_ens = df_ens__[df_ens__['Gene name'] == gene.upper()]
            if not df_ens.empty:
                mir_start = df_mir.loc[idx, 'start']
                mir_end = df_mir.loc[idx, 'end']
                chromosome = df_mir.loc[idx, 'chromosome']

                ens_tss = df_ens.iloc[0]['Transcription start site (TSS)']
                strand = df_ens.iloc[0]['Strand']
                ens_chromosome = 'chr' + df_ens.iloc[0]['Chromosome/scaffold name']

                if chromosome != ens_chromosome:
                    print(gene, chromosome, ens_chromosome)

                if strand > 0:
                    mir_tss = mir_start
                else:
                    mir_tss = mir_end
                contents.append([gene, abs(mir_tss - ens_tss)])

        df_rep = pd.DataFrame(contents, columns=['mirna', 'distance'])
        df_rep.to_excel('compare_fantom_ensembl.xlsx', index=None)

    def compare_mir(self):
        fpath = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export.db')
        con = sqlite3.connect(fpath)
        df_ens__ = pd.read_sql_query("SELECT * FROM 'Ensembl'", con)

        contents = []
        df_mir = pd.read_excel('miRNA.xlsx')
        for idx in df_mir.index:
            print('{} / {}'.format(idx + 1, df_mir.shape[0]))
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
                ens_chromosome = 'chr' + df_ens.iloc[0]['Chromosome/scaffold name']

                if chromosome != ens_chromosome:
                    print(mirna, chromosome, ens_chromosome)

                if strand > 0:
                    mir_tss = mir_start
                else:
                    mir_tss = mir_end
                contents.append([mirna, abs(mir_tss - ens_tss)])

        df_rep = pd.DataFrame(contents, columns=['mirna', 'distance'])
        df_rep.to_excel('compare_fantom_ensembl.xlsx', index=None)

    def match(self):
        fpath = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export.db')
        con = sqlite3.connect(fpath)
        df_38 = pd.read_sql_query("SELECT * FROM 'Ensembl'", con)

        fpath = os.path.join(self.root, 'database/ensembl/TSS', 'Ensembl_hg19.filter')
        df_19 = pd.read_csv(fpath, sep='\t', names=['chr_old', 'start_old', 'end_old', 'arrow', 'chr_new', 'start_new', 'end_new'])
        df_19.index = df_19['chr_old'].astype(str) + ';' + df_19['start_old'].astype(str) + ';' + df_19['end_old'].astype(str)
        df_19 = df_19[~df_19.index.duplicated(keep='first')]

        drop_list = []
        for idx in df_38.index:
            if idx % 1000 == 0 or idx == df_38.shape[0] - 1:
                print('{:0.2f}%'.format(100 * (idx + 1) / df_38.shape[0]))
            start_old = df_38.loc[idx, 'Transcript start (bp)']
            end_old = df_38.loc[idx, 'Transcript end (bp)']
            chromosome_old = 'chr' + df_38.loc[idx, 'Chromosome/scaffold name']
            key = ';'.join([str(chromosome_old), str(start_old), str(end_old)])

            if key not in df_19.index:
                drop_list.append(idx)
                continue
            chromosome_new = df_19.loc[key, 'chr_new']
            start_new = df_19.loc[key, 'start_new']
            end_new = df_19.loc[key, 'end_new']

            df_38.loc[idx, 'Transcript start (bp)'] = start_new
            df_38.loc[idx, 'Transcript end (bp)'] = end_new
            df_38.loc[idx, 'Chromosome/scaffold name'] = chromosome_new

        df_38 = df_38.drop(drop_list)
        out_path = os.path.join(self.root, 'database/ensembl/TSS', 'mart_export_hg19.csv')
        df_38.to_csv(out_path, sep='\t', index=None)


if __name__ == '__main__':
    ee = Extract_Ensembl()
    ee.match()
