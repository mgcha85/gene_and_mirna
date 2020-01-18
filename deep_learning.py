import pandas as pd
import sqlite3
import os
import socket


class deep_learning:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def model(self, X, Y):
        from keras.models import Sequential, load_model
        from keras.layers import Dense, Dropout, Activation
        from keras.optimizers import SGD

        model = Sequential()
        model.add(Dense(20, activation="tanh", input_dim=5, kernel_initializer='uniform'))
        model.add(Dense(1, activation='linear', kernel_initializer='uniform'))

        model.compile(loss='mse', optimizer='adam', metrics=['accuracy'])
        model.fit(X, Y, epochs=100, batch_size=10, verbose=2)

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/cell_lines/out', 'regression_100_nz.db')
        con = sqlite3.connect(fpath)

        df_x = pd.read_sql("SELECT * FROM 'X'", con, index_col='tid')
        df_y = pd.read_sql("SELECT * FROM 'Y'", con, index_col='miRNA')

        self.model(df_x, df_y)


if __name__ == '__main__':
    dl = deep_learning()
    dl.run()