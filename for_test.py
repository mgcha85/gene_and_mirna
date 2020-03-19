import sqlite3
import pandas as pd
import os
from Database import Database

fpath = os.path.join('D:/Bioinformatics/database/Fantom/v5/tissues/CAGE_tag_tissue_spt.db')
con = sqlite3.connect(fpath)
tlist = Database.load_tableList(con)
tissues = set([x.split('_')[0] for x in tlist])

# for tis in sorted(list(tissues)):
#     df = pd.read_sql("SELECT score FROM '{}_chr9_+' WHERE start<107510064 AND end>107509864".format(tis), con)
#     print(tis)
#     print(df['score'])


fpath = os.path.join('D:/Bioinformatics/database/RNA-seq/out/RNA_seq_tissue.db')
con = sqlite3.connect(fpath)
tlist = Database.load_tableList(con)
tissues = set.intersection(set(tlist), tissues)

for tis in sorted(list(tissues)):
    df = pd.read_sql("SELECT FPKM FROM '{}' WHERE reference_id='ENST00000374767.5_2'".format(tis), con)
    print(tis)
    print(df['FPKM'])
