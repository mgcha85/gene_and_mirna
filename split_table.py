import os
import pandas as pd
import socket
import sqlite3


class Split_table:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def categorizing(self):
        with open('cell_lines.txt', 'rt') as f:
            contents = f.read()
        contents = contents.split('\n')
        data = [[x.split('%20')[0][7:], x] for x in contents]
        df = pd.DataFrame(data=data, columns=['key', 'columns'])
        df_grp = df.groupby('key')
        return df_grp

    def run(self):
        df_grp = self.categorizing()
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc_full.db')
        con = sqlite3.connect(fpath)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df = pd.read_csv(fpath)
        for group, df_cols in df_grp:
            if len(group) == 0:
                continue

            try:
                columns = ['chromosome', 'start', 'end', 'strand'] + list(df_cols['columns'].values)
                df_cell = df[columns]
                df_cell = df_cell[df_cell.iloc[:, 3:].sum(axis=1) > 0]
                df_cell[columns].to_sql(group, con, if_exists='replace', index=None)
            except Exception as e:
                print(e)
                continue


if __name__ == '__main__':
    st = Split_table()
    st.run()
