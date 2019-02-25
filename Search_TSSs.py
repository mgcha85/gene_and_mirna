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

    def split_table_by_cell_lines(self, cell_lines):
        fname__ = 'hg19.cage_peak_phase1and2combined_counts.osc'
        dirname = os.path.join(self.root, 'database/Fantom/v5')
        fpath = os.path.join(dirname, fname__ + '.csv')
        df = pd.read_csv(fpath)
        # columns = list(df.columns)
        # with open('columns.txt', 'wt') as f:
        #     f.write('\n'.join(columns))
        # exit(1)

        for cline in cell_lines:
            print(cline)
            columns = ['chromosome', 'start', 'end', 'strand']
            for col in df.columns:
                if cline in col:
                    columns.append(col)

            out_path = os.path.join(dirname, fname__ + '_{}.csv'.format(cline))
            df[columns].to_csv(out_path, index=None)


if __name__ == '__main__':
    st = Search_TSSs()
    cell_lines = ['K562', 'HepG2', 'A549', 'GM12878', 'HEK293', 'MCF-7']
    st.split_table_by_cell_lines(cell_lines)
