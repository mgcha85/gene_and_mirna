import sqlite3
import os
import socket
import pandas as pd
import numpy as np


class comparison:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.table_names = {}

    def run(self):
        import matplotlib.pyplot as plt

        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz_others.db')
        con = sqlite3.connect(fpath)

        df_cor = pd.read_sql("SELECT * FROM 'corr'", con, index_col='transcript_id')
        df_coef = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='tid')

        nz = np.zeros(df_coef.shape[1])
        for i, col in enumerate(df_coef.columns):
            nz[i] = len(df_coef[col][df_coef[col] != 0])

        plt.subplot(111)
        xaxis = range(df_cor.shape[1])
        plt.plot(xaxis, df_cor.mean())
        plt.xticks(xaxis[::10], df_cor.columns[::10], fontsize=6, rotation=30)
        plt.title('corr coeff by a miRNA')
        # plt.subplot(212)
        # plt.plot(range(df_coef.shape[1]), nz / df_coef.shape[0])
        plt.grid()
        plt.show()


if __name__ == '__main__':
    comp = comparison()
    comp.run()
