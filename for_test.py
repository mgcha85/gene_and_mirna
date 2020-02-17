import sqlite3
import pandas as pd

fpath = '/home/mingyu/Bioinformatics/database/target_genes/predictions_processed.db'
con = sqlite3.connect(fpath)
df_mir = pd.read_sql("SELECT * FROM 'miranda_grp_mir'", con, index_col='mir')
df_rna = pd.read_sql("SELECT * FROM 'rna22_grp_mir'", con, index_col='mir')
df_ts = pd.read_sql("SELECT * FROM 'ts_grp_mir'", con, index_col='mir')

genes1 = df_mir.loc['hsa-let-7a-3p', 'genes'].split(';')
genes2 = df_rna.loc['hsa-let-7a-3p', 'genes'].split(';')
genes3 = df_ts.loc['hsa-let-7a-3p', 'genes'].split(';')

print('see')
