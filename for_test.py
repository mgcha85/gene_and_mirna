import pandas as pd
import os
import sqlite3
from Database import Database

con = sqlite3.connect('/home/mingyu/Bioinformatics/database/consistent_miRNA_330.db')
con_out = sqlite3.connect('/home/mingyu/Bioinformatics/database/consistent_miRNA_330_2.db')
tlist = Database.load_tableList(con)
for tname in tlist:
    df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='miRNA')
    df = df.loc[~df.index.duplicated(keep='first')]
    df.to_sql(tname, con_out, if_exists='replace')
