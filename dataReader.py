#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:13:45 2021

@author: BrianTCook
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')

from galpy.util import bovy_conversion
from galpy.potential import MWPotential2014

np.random.seed(42)
    
def jacobiRadius(m, r, z):
    
    galMass = 0.
    
    for massElement in MWPotential2014:

        galMass += massElement.mass(r, z= z * bovy_conversion.mass_in_msol(220., 8.)) #MSun

    return galMass

def streamPlotter(file):
    
    df = pd.read_csv(file, sep=' ', names=['mass','x','y','z','vx','vy','vz'])
    dfStreams = np.array_split(df, 3)

    fileStrs = file.split('_')
    j = fileStrs[2]

    fig, axs = plt.subplots(1, 3, figsize=(8, 3), constrained_layout=True)
    
    for i, stream in enumerate(dfStreams):
        
        axi = axs[i]
        
        yVals = [ y/1000. for y in stream['x'].tolist() ]
        zVals = [ z/1000. for z in stream['y'].tolist() ]
        axi.scatter(yVals, zVals, c='k', s=2)
        axi.scatter(yVals[0], zVals[0], c='r', s=12)
        
    fig.text(0.5, -0.05, r'$y_{\rm MW}$ (kpc)', ha='center', fontsize=16)
    fig.text(-0.05, 0.6, r'$z_{\rm MW}$ (kpc)', rotation='vertical', fontsize=16)
        
    plt.tight_layout()
    plt.savefig('snapshot_%s.pdf'%(j), bbox_inches = "tight")
    
if __name__ in '__main__':
    
    print('hello world')
    files = glob.glob('*.ascii')
    streamPlotter(files[16])