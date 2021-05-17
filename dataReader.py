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

import time

#plt.rc('text', usetex = True)
#plt.rc('font', family = 'serif')

from galpy.util import bovy_conversion
from galpy.potential import MWPotential2014

np.random.seed(42)
    
def jacobiRadius(m, r, z):
    
    galMass = 0.

    for massElement in MWPotential2014:
        galMass += massElement.mass(r, z=z) * bovy_conversion.mass_in_msol(220., 8.) #MSun
        
    galactocentricDistance = np.sqrt(r**2. + z**2.)

    return galactocentricDistance * (m / galMass)**(1/3.)

def tidalRadii(nStreams, nSnapshots, files, streamNames):
    
    t0 = time.time()
    outputArr = np.zeros((nSnapshots, nStreams))
    tVals = np.linspace(0., 35., 36)
    
    for i, file in enumerate(files):
        
        print('')
        print('currently at snapshot %.0f Myr: '%(tVals[i]))
        print('wall time: %.02f minutes'%((time.time()-t0)/60.))
        
        df = pd.read_csv(file, sep=' ', names=['mass','x','y','z','vx','vy','vz'])
        dfBigParticles = df[df['mass'] > 1.].reset_index()
        #dfBigParticles = dfBigParticles.reset_index()
        
        for j, row in dfBigParticles.iterrows():
            
            xKPC, yKPC, zKPC = row['x']/1000., row['y']/1000., row['z']/1000.
            particleMass = row['mass']
            
            particleR, particleZ = np.sqrt(xKPC**2. + yKPC**2.), zKPC
    
            outputArr[i,j] = jacobiRadius(particleMass, particleR, particleZ) * 1000.
    
    for k in range(nStreams):
        plt.semilogy(tVals, outputArr[:,k])
    
    plt.annotate(r'R_{\rm max}$', xy = (0.05, 0.95), xycoords = 'axes fraction', fontsize = 8)
    
    return outputArr

'''
def streamPlotter(nStreams, nSnapshots, files):
    
    ax1.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
    ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
    ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
    ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
    ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

    for file in files:

		df = pd.read_csv(file, sep=' ', names=['mass','x','y','z','vx','vy','vz'])
		dfStreams = np.array_split(df, nStreams)

		fileStrs = file.split('_')
		j = fileStrs[2]

		fig, axs = plt.subplots(2, 2, constrained_layout=True)

		for i, stream in enumerate(dfStreams):

			yVals = [ y/1000. for y in stream['x'].tolist() ]
			zVals = [ z/1000. for z in stream['y'].tolist() ]
			axi.scatter(yVals, zVals, c='k', s=2)
			axi.scatter(yVals[0], zVals[0], c='r', s=12)

		fig.text(0.5, -0.05, r'$y_{\rm MW}$ (kpc)', ha='center', fontsize=16)
		fig.text(-0.05, 0.6, r'$z_{\rm MW}$ (kpc)', rotation='vertical', fontsize=16)

		plt.tight_layout()
		plt.savefig('snapshot_%s.pdf'%(j), bbox_inches = "tight")
		plt.close()
'''
    
if __name__ in '__main__':
    
    nStreams, nSnapshots = 5, 36
    
    files = glob.glob('*.ascii')
    tidalRadii(nStreams, nSnapshots, files, streamNames)
    
    #streamPlotter(files, nStreams)
