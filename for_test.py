from Database import Database
import os
import sqlite3
import pandas as pd

df = pd.read_excel('/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/Papers/complement/supp_gkv608_nar-00656-h-2015-File009.xls', header=1)
print('well')
root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
dirname = os.path.join(root, 'database/Fantom/v5/tissues/out')
flist = os.listdir(dirname)

out_path = os.path.join(root, 'database/Fantom/v5/tissues', 'hg19.final_confirmed_tss.db')
out_con = sqlite3.connect(out_path)

for fname in flist:
    fpath = os.path.join(dirname, fname)
    con = sqlite3.connect(fpath)
    tlist = Database.load_tableList(con)

    for tname in tlist:
        df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
        df.to_sql(tname, out_con, if_exists='replace', index=None)

