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

df = pd.read_excel('correlation_500.xlsx')

for idx in df.index:
    scores = df.loc[idx, 'scores']
    scores = scores.split(':')
    print(scores)
