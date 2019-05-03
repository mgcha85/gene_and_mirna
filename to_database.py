import sqlite3
import pandas as pd
import socket
import os


class To_database:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    def gtf_to_db(self):
        fname = 'gencode.v30lift37.annotation'
        fpath = os.path.join(self.root, 'database/gencode', '{}.db'.format(fname))
        out_con = sqlite3.connect(fpath)

        chunksize = 1 << 29     # 512 MB
        for df_chunk in pd.read_csv(fpath.replace('.db', '.gtf'), sep='\t', chunksize=chunksize, names=self.gtf_columns, comment='#'):
            df_chunk.to_sql(os.path.splitext(fname)[0], out_con, if_exists='append', index=None)


if __name__ == '__main__':
    td = To_database()
    td.gtf_to_db()

