import os
import pandas as pd
import sqlite3
from Database import Database
import gzip


dirname = "D:/Bioinformatics/database/Fantom/v5/cell_lines"
fpath = os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz')
con = sqlite3.connect(os.path.join(dirname, 'hg19.cage_peak_phase1and2combined_tpm_ann.osc.db'))
for df_sub, comments in pd.read_csv(fpath, sep='\t', iterator=True, chunksize=1 << 10, comment='#', return_comments=True):
    df_sub.to_excel('temp.xlsx', index=None)
    with open('temp.txt', 'wt') as f:
        f.write(comments)
    break

  compressor = gzip.GzipFile(fileobj=stream, mode='w')

with gzip.open(fpath, 'rt') as f:
    for line in f:
        print('got line', line)

#
# for fname in ["sum_fan_gene_100.db", "sum_fan_mir_100.db"]:
#     fpath = os.path.join(dirname, fname)
#     con = sqlite3.connect(fpath)
#     dfs = []
#     for tname in Database.load_tableList(con):
#         dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
#     df = pd.concat(dfs)
#     df.to_excel(fpath.replace('.db', '.xlsx'), index=False)
