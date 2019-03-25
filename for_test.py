import pandas as pd
import sqlite3
from Database import Database


fpath = '/home/mingyu/Dropbox/Tss_map_table.db'
con = sqlite3.connect(fpath)
tlist = Database.load_tableList(con)

out_path = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/Papers/Tss_map/Tss_map_table2.db'
out_con = sqlite3.connect(out_path)
for tname in tlist:
    df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
    df.to_sql(tname, out_con, if_exists='replace', index=None)
