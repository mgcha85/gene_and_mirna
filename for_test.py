# import numpy as np
# from scipy.stats import wilcoxon
#
# x = [0.00054926, 0.00101568, 0.00088054, 0.0027738, 0.00096129, 0.00217144, 0.00099867, 0.00594059, 0.00135477, 0.00107772, 0.00333265, 0.00110879, 0.00228319, 0.00470631, 0.00115906, 0.00504211, 0.00767846,0.00074081, 0.00074581, 0.00438994, 0.00149771, 0.0060626, 0.00072757, 0.00674521]
# y = [1.63984243e-03, 7.40561331e-05, 1.23005911e-03, 2.29219876e-03, 3.14356099e-03, 1.35880550e-03, 4.90689694e-03, 5.95341709e-03, 1.33250515e-04, 7.63679485e-04, 4.88740198e-03, 1.80198209e-03, 1.58698933e-03, 4.96086045e-03, 4.58201903e-03, 4.74290841e-03, 3.89539991e-03, 1.45752779e-03, 2.77170225e-03, 3.95982779e-03, 1.89376775e-04, 3.48414076e-03, 2.04057036e-03, 4.27119215e-03]
#
# x = np.array(x)
# y = np.array(y)
#
# print((x - y).sum())
# T, p_value = wilcoxon(x, y, zero_method="wilcox")
# print('W = {:0.0f}, p-value = {:0.4f}'.format(T, 1 - p_value / 2))


import sqlite3
import pandas as pd
import numpy as np
import os
from Database import Database

con = sqlite3.connect("D:/Bioinformatics/database/Fantom/v5/tissues/sum_fan_gene_100.db")
dfs = []
for tname in Database.load_tableList(con):
    df = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
    dfs.append(df.shape[0])

print(sum(dfs))

# dirname = 'D:/Bioinformatics/database/Fantom/v5/tissues'
# df = pd.read_csv(os.path.join(dirname, '00_human.tissue.hCAGE.hg19.assay_sdrf.txt'), sep='\t')
# df['Comment [sample_name]'] = df['Comment [sample_name]'].str.split(' - ').str[0]
# df['Comment [sample_name]'] = df['Comment [sample_name]'].str.split(', ').str[0]
#
# contents = []
# for grp, df_sub in df.groupby('Comment [sample_name]'):
#     contents.append([grp, ';'.join(df_sub['File Name.1']), df_sub.shape[0]])
# df_res = pd.DataFrame(contents, columns=['tissue', 'file names', '#'])
# df_res.to_excel('file_list.xlsx', index=None)

# df = pd.read_excel('file_list.xlsx')
# con = sqlite3.connect(os.path.join(dirname, 'FANTOM_tissue.db'))
# tlist = Database.load_tableList(con)
# remain = set(df['tissue']) - set(tlist)
# print(remain)