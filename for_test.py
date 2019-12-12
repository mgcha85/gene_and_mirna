import pandas as pd
import sqlite3
from Database import Database

con = sqlite3.connect("/home/mingyu/Bioinformatics/database/Fantom/v5/cell_lines/sum_fan_gene_100.db")
cnt = 0
for tname in Database.load_tableList(con):
    df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='transcript_id')
    cnt += df.shape[0]
print(cnt)
