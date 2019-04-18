import pandas as pd
import numpy as np
from scipy.stats import spearmanr

n_rows = 2500
cols = ['a', 'b', 'c', 'd', 'e', 'f', 'g']

df = pd.DataFrame(np.random.random(size=(n_rows, len(cols))), columns=cols)
v = np.random.random(size=len(cols))

# original implementation
corr, _ = zip(*df.apply(lambda x: spearmanr(x,v), axis=1))
corr = pd.Series(corr)

# modified implementation
df1 = df.rank(axis=1)
v1 = pd.Series(v, index=df.columns).rank()
corr1 = df1.corrwith(v1, axis=1)

print(corr1)