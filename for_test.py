import sqlite3
import pandas as pd


fpath = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/Fantom/v5/hg19.cage_peak_phase1and2combined_counts.osc.csv'
df = pd.read_csv(fpath)

df = df[['chromosome', 'start', 'end', 'strand']]
con = sqlite3.connect(fpath.replace('.csv', '.db'))
df.to_sql('Fantom_count', con, if_exists='replace')
