import os
import pandas as pd
import sqlite3
from Database import Database

dirname = "D:/Bioinformatics/database/Fantom/v5/cell_lines"

for fname in ["sum_fan_gene_100.db", "sum_fan_mir_100.db"]:
    fpath = os.path.join(dirname, fname)
    con = sqlite3.connect(fpath)
    dfs = []
    for tname in Database.load_tableList(con):
        dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
    df = pd.concat(dfs)
    df.to_excel(fpath.replace('.db', '.xlsx'), index=False)
