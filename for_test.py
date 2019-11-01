# import pandas as pd
# import os
# import sqlite3
# from Database import Database
#
# con = sqlite3.connect('/home/mingyu/Bioinformatics/database/consistent_miRNA_330.db')
# con_out = sqlite3.connect('/home/mingyu/Bioinformatics/database/consistent_miRNA_330_2.db')
# tlist = Database.load_tableList(con)
# for tname in tlist:
#     df = pd.read_sql("SELECT * FROM '{}'".format(tname), con, index_col='miRNA')
#     df = df.loc[~df.index.duplicated(keep='first')]
#     df.to_sql(tname, con_out, if_exists='replace')

starts = [75062507,75062545,75062636,75062663,75062704,75062773,75062812,75062873]
ends = [75062507,75062545,75062789,75062663,75062704,75062973,75062968,75062874]

smin = max(starts)
emin = max(ends)
print(max([smin, emin]))
print(int((75062507 + 75062973) // 2))