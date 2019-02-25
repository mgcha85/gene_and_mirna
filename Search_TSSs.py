import pandas as pd
import socket
import os


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

    def load_cnt_table(self):
        fpath = os.path.join(self.root, 'database', 'Fantom', 'v5', 'hg19.cage_peak_phase1and2combined_counts.osc.csv')
        df = pd.read_csv(fpath)
        columns = list(df.columns)
        contents = '\n'.join(columns)
        with open('cell_lines.txt', 'wt') as f:
            f.write(contents)


if __name__ == '__main__':
    st = Search_TSSs()
    st.load_cnt_table()
