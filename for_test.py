import pandas as pd

df = pd.read_excel('/home/mingyu/Documents/correlation.xlsx')
df.T.to_excel('/home/mingyu/Documents/correlation2.xlsx')