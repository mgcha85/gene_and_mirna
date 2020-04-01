import pandas as pd
import sqlite3
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

dirname = "D:/Users/mgcha/Downloads"
# fpaths = [os.path.join(dirname, "41467_2018_7746_MOESM7_ESM.txt")]

con = sqlite3.connect(os.path.join(dirname, "41467_2018_7746_MOESM.db"))
# for fpath in fpaths:
#     df = pd.read_csv(fpath, sep='\t')
#     loc = df['SNP(hg19)'].str.split(':', expand=True)
#     loc.columns = ['chromosome', 'start']
#     df = pd.concat([df, loc], axis=1)
#     df = df.drop(['SNP(hg19)'], axis=1)
#     df['start'] = df['start'].astype(int)
#     df['end'] = df['start'] + 1
#     columns = list(df.columns)
#     score = columns.pop(0)
#     columns.append(score)
#     df = df[columns]
#     dirname, fname = os.path.split(fpath)
#     df.to_sql(os.path.splitext(fname)[0], con, index=None, if_exists='replace')
# exit(1)

chromosome = 'chr9'
start = 127047500
end = start + 2000

titles = ['active regions', 'tiled regions', 'high-resolution driver elements',
          'SNPs with allele-specific activity (raw p-value < 0.05)']
n = 4
_, ax = plt.subplots(n)
for i in range(n):
    tname = '41467_2018_7746_MOESM{}_ESM'.format(i+4)
    df = pd.read_sql("SELECT * FROM '{}' WHERE chromosome='{}' AND start<{} AND end>{}".format(tname, chromosome, end, start), con)
    for idx in df.index:
        buffer = np.zeros(end - start)
        sidx = df.loc[idx, 'start'] - start
        eidx = df.loc[idx, 'end'] - df.loc[idx, 'start'] + sidx

        buffer[sidx:eidx] = 1
        xaxis = np.arange(buffer.shape[0])
        ax[i].plot(xaxis, buffer)

        # ax[i].set_title(titles[i])
        ax[i].set_xticks(xaxis[::400])
        ax[i].set_xticklabels(['{}:{}'.format(chromosome, x) for x in range(start, end, 400)], fontsize=6)
plt.show()