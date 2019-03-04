import sqlite3
from Database import Database
import pandas as pd


fpath = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/Papers/Tss_map/Tss_map.db'
con = sqlite3.connect(fpath)
tlist = Database.load_tableList(con)
cell_line = 'HEK293'

out_fpath = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/Papers/Tss_map/Tss_map_table2.db'
out_con = sqlite3.connect(out_fpath)
tlist = [x for x in tlist if cell_line in x]
for tname in tlist:
    df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
    df.to_sql(tname, out_con, if_exists='replace', index=None)
