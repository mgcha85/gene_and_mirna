import pandas as pd
import os
import socket


class Lasso:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinfomatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def trim_genes(self, genes):
        result = []
        for gene in genes:
            if '///' in gene:
                result.append(gene.split('///')[0])
            else:
                result.append(gene)
        return result

    def load_import_genes(self):
        df = pd.read_csv(os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'important_genes.txt'), sep='\t', names=['miRNA', 'GENEs', 'Weights'], index_col=0)
        df = df.dropna(subset=['GENEs'], axis=0)
        for mir in df.index:
            gene_set = df.loc[mir, 'GENEs'].split(';')
            gene_set = self.trim_genes(gene_set)
            gene_set = ';'.join(gene_set)
            df.loc[mir, 'GENEs'] = gene_set
        df.to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'important_genes.xlsx'))

    def set_import_genes(self):
        contents = []
        df = pd.read_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'important_genes.xlsx'), index_col=0)
        for mir in df.index:
            gene_set = df.loc[mir, 'GENEs'].split(';')
            weight_set = df.loc[mir, 'Weights'].split(';')
            for g, w in zip(gene_set, weight_set):
                contents.append([mir, g, w])
        df_res = pd.DataFrame(contents, columns=['miRNA', 'gene', 'weight'])
        df_res.to_excel(os.path.join(self.root, 'database/Fantom/v5/cell_lines', 'important_genes2.xlsx'))


if __name__ == '__main__':
    lasso = Lasso()
    lasso.set_import_genes()
