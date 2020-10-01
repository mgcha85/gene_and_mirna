import sqlite3
import pandas as pd
from Database import Database

con = sqlite3.connect("D:/Bioinformatics/database/gencode/high_correlated_fan_rna_100.db")
cnt = 0
for tname in Database.load_tableList(con):
	df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='transcript_id')
	cnt += df.shape[0]
print(cnt)
