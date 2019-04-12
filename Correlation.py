import socket
import os
import sqlite3
import pandas as pd


class Correlation:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def load_fantom(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.txt')
        df = pd.read_csv(fpath, comment='#')
        return df
        
    def run(self):
        con = sqlite3.connect(os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.db'))
        df_fantom = self.load_fantom()
        df_fantom.to_sql('hg19.cage_peak_phase1and2combined_counts.osc', con, if_exists='replace', index=None)

        # fpath = os.path.join(self.root, 'database', 'genecode.db')
        # con = sqlite3.connect(fpath)


if __name__ == '__main__':
    cor = Correlation()
    cor.run()
