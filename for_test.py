import pandas as pd


fpath = '/lustre/fs0/home/mcha/Bioinformatics/database/Fantom/v5/hg19.cage_peak_phase1and2combined_counts.osc.csv'
df = pd.read_csv(fpath)
df = df[['chromosome', 'start', 'end', 'strand']]
df.to_csv(fpath.replace('.csv', '_loc.csv'), index=None)
