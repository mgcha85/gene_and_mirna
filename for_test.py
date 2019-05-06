import sqlite3
import pandas as pd
from Database import Database


fpath = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/ensembl/TSS/mart_export_hg19.db'
con = sqlite3.connect(fpath)
out_con = sqlite3.connect(fpath.replace('.db', '_pc.db'))

tlist = Database.load_tableList(con)
for tname in tlist:
    _, chr, str = tname.split('_')

    df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
    df = df[df['Transcript type'] == 'protein_coding']
    if str == '+':
        df['Transcription start site (TSS)'] = df['start']
    else:
        df['Transcription start site (TSS)'] = df['end']
    df.to_sql(tname, out_con, index=None, if_exists='replace')
