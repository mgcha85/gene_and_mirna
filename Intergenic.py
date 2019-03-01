import sqlite3
import os
import socket
import pandas as pd


class Intergenic:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Independent_research/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        elif 'evc' in hostname:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        else:
            print('wrong option')
            return

    def run(self):
        dirname = os.path.join(self.root, 'database')
        fname = 'fantom5.db'

        fpath = os.path.join(dirname, fname)
        con = sqlite3.connect(fpath)

        df = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates' WHERE Type='intronic'", con)

        for idx in df.index:
            promoter = df.loc[idx, 'promoter']
            if promoter.split('@')[0] == 'p':
                df.loc[idx, 'type'] = 'intergenic'
            else:
                df.loc[idx, 'type'] = 'intronic'

        df.to_sql('human_promoters_wo_duplicates', con, if_exists='replace', index=None)


if __name__ == '__main__':
    int = Intergenic()
    int.run()
