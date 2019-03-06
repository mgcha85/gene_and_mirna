import os
import pandas as pd
import socket
import sqlite3


class Split_table:
    def __init__(self):
        hostname = socket.gethostname()
        print(hostname)
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
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
        keys = set([x.split('%20')[0][7:] for x in contents])

        map = {}
        for k in keys:
            if len(k) == 0:
                continue
            for con in contents:
                if k in con:
                    if k in map:
                        map[k].append(con)
                    else:
                        map[k] = [con]
        return map

    def run(self):
        map = self.categorizing()
        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19_cage_peak_phase1and2combined_counts_osc.db')
        con = sqlite3.connect(fpath)

        fpath = os.path.join(self.root, 'database/Fantom/v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df = pd.read_csv(fpath)
        for key, cols in map.items():
            df[cols].to_sql(key, con, if_exists='replace', index=None)


if __name__ == '__main__':
    st = Split_table()
    st.run()