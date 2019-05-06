import sqlite3
import pandas as pd


fpath = 'D:/Bioinformatics/database/Fantom/v5/hg19.cage_peak_phase1and2combined_counts.osc.csv'
df = pd.read_csv(fpath)

df = df[['chromosome', 'start', 'end', 'strand']]
con = sqlite3.connect(fpath.replace('.csv', '.db'))
df.to_sql('Fantom', con, if_exists='replace', index=None)
