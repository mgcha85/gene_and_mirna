import sqlite3
import pandas as pd
import socket
import os


class To_database:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    def gtf_to_db(self):
        fname = 'genes'
        fpath = os.path.join(self.root, 'software/stringtie-1.3.3b', '{}.db'.format(fname))
        out_con = sqlite3.connect(fpath)

        chunksize = 1 << 29     # 512 MB
        for df_chunk in pd.read_csv(fpath.replace('.db', '.gtf'), sep='\t', chunksize=chunksize, names=self.gtf_columns, comment='#'):
            df_chunk['chromosome'] = 'chr' + df_chunk['chromosome'].astype('str')
            df_chunk_chr = df_chunk.groupby('chromosome')
            for chr, df_chunk_sub in df_chunk_chr:
                df_chunk_chr_str = df_chunk_sub.groupby('strand')
                for str, df_chunk_sub_sub in df_chunk_chr_str:
                    df_chunk_sub_sub.to_sql(os.path.splitext(fname)[0] + '_{}_{}'.format(chr, str), out_con, if_exists='append', index=None)

    def temp(self):
        fpath = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/Papers/mirSTP/GSM1006730_Gro-Seq_1hr.bed'
        out_path = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/Papers/mirSTP/GSM1006730_Gro-Seq_1hr_out.bed'
        columns = ['chromomsome', 'start', 'end', 'attribute', 'score']

        df = pd.read_csv(fpath, sep='\t', names=columns, comment='#')
        df.loc[:, 'strand'] = '.'
        df.to_csv(out_path, sep='\t', index=None, header=False)


if __name__ == '__main__':
    td = To_database()
    td.temp()

