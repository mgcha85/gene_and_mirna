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

fpath1 = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/histone_features.db'
fpath2 = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/histone_features2.db'

con1 = sqlite3.connect(fpath1)
con2 = sqlite3.connect(fpath2)

tlist = Database.load_tableList(con2)

for tname in tlist:
    df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con2)
    df.to_sql(tname, con1, index=None)
