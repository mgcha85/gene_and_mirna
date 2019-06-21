import pandas as pd
import sqlite3
import os
import socket
from Database import Database
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Server import Server
import sys
import pickle as pkl
import numpy as np


# df = pd.read_csv('/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/miRNA_redefined_Tss_coor_4.bed', sep='\t', index_col=0)
# con = sqlite3.connect('/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/fantom5.db')
# df_fan = pd.read_sql_query("SELECT * FROM 'human_promoters_wo_duplicates'", con, index_col='premiRNA')
# mir = set.intersection(set(df.index), set(df_fan.index))
# df_fan.loc[mir].to_sql('validate_mir_tss', con, if_exists='replace')
# exit(1)

df_100 = pd.read_excel(
    '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/Fantom/v5/tissues/out/fantom_cage_by_tissue_100.xlsx',
    index_col=0)
df_500 = pd.read_excel(
    '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/Fantom/v5/tissues/out/fantom_cage_by_tissue_500.xlsx',
    index_col=0)

genes = list(set(df_100.index) - set(df_500.index))

with open('why.txt', 'wt') as f:
    f.write('\n'.join(sorted(genes)))
exit(1)

fpath = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/fantom5.db'
con = sqlite3.connect(fpath)
out_con = sqlite3.connect(fpath.replace('.db', '_mir.db'))
df = pd.read_sql_query("SELECT * FROM 'validate_mir_tss'", con)
df_chr = df.groupby('chromosome')
for chr, df_sub in df_chr:
    df_str = df_sub.groupby('strand')
    for str, df_sub_sub in df_str:
        df_sub_sub = df_sub_sub.rename(columns={'premiRNA': 'gene_name'})
        df_sub_sub['start'] = df_sub_sub['tss']
        df_sub_sub['end'] = df_sub_sub['tss']
        df_sub_sub.drop('tss', axis=1).to_sql('human_promoters_{}_{}'.format(chr, str), out_con, if_exists='replace', index=None)
exit(1)


dirname = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/Fantom/v5/tissues'
flist = os.listdir(dirname)
tissues = []
for fname in flist:
    if os.path.isdir(os.path.join(dirname, fname)) is True and fname != 'out':
        tissues.append(fname)

df = pd.Series(data=tissues)
df.to_csv('FANTOM_tissue.csv', index=None)
