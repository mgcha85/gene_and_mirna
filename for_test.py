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

con = sqlite3.connect('/home/mingyu/Bioinformatics/database/Fantom/v5/cell_lines/out/cross_validation/nz/regression_100_0_test.db')
df_x = pd.read_sql("SELECT * FROM 'X'", con, index_col='tid').T
df_y = pd.read_sql("SELECT * FROM 'Y'", con, index_col='miRNA').T
df_b = pd.read_sql("SELECT * FROM 'coefficient'", con, index_col='tid').T

df_x.to_csv('/home/mingyu/Bioinformatics/database/Fantom/v5/cell_lines/out/cross_validation/nz/X.csv', header=False, index=None, sep='\t')
df_y.to_csv('/home/mingyu/Bioinformatics/database/Fantom/v5/cell_lines/out/cross_validation/nz/Y.csv', header=False, index=None, sep='\t')
df_b.to_csv('/home/mingyu/Bioinformatics/database/Fantom/v5/cell_lines/out/cross_validation/nz/B.csv', header=False, index=None, sep='\t')
