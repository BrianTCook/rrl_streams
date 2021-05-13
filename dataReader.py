#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:13:45 2021

@author: BrianTCook
"""

import pandas as pd
import matplotlib.pyplot as plt
import glob

files = glob.glob('*.ascii')

print(files)

for file in files:

    df = pd.read_csv(file, sep=' ', names=['mass','x','y','z','vx','vy','vz'])
    plt.scatter(df['x']/1000.,df['y']/1000.)
    plt.xlim(-50., 50.)
    plt.ylim(-50., 50.)
    plt.show()