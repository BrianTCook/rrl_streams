#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 09:21:07 2021

@author: BrianTCook
"""

import numpy as np
import pandas as pd
import time
import math
import glob
import matplotlib.pyplot as plt

plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')

from astropy.coordinates import SkyCoord
import astropy.units as u
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
from galpy.potential import MWPotential2014, to_amuse

np.random.seed(42)

projectDir = '/Users/BrianTCook/Desktop/MITLL_project/clustering_rr_lyrae/'

#known globular clusters
columns_gc = ['ID', 'Name', 'RAh', 'RAm', 'RAs', 'DECd', 'DECm', 'DECs', 'L', 'B', 'R_Sun',  'R_gc', 'X', 'Y', 'Z']
dfGCs = pd.read_csv(projectDir+'milkyway_gcs_data.csv', sep=',', names=columns_gc)
dfStars = pd.read_csv(projectDir+'RRLs_with_Gaia.csv', names=[  'dRR', 'ddRR', 'RA', 'DEC', 'l', 'b', 
                                                              'mV', 'dmV', 'MV', 'dMV', 'AV', 'dAV', 
                                                              'xcoord', 'ycoord', 'zcoord', 'source', 'muRA', 'muDEC' ]) #np.loadtxt('X_AAVSO.txt')
dfGroupings = pd.read_csv(projectDir+'forest_groupings.csv')
dfNewGCs = dfGroupings[dfGroupings['Cluster Tag']=='GC']
dfOthers = dfGroupings[dfGroupings['Cluster Tag']=='Other']
dfSSs = dfGroupings[dfGroupings['Cluster Tag']=='SS']

specialGCs = [ 'NGC 5904', 'NGC 5139', 'NGC 6121' ]
gcMasses = { 'NGC 5904':3.72e5, 'NGC 5139':3.55e6, 'NGC 6121':9.69e4}

dfSpecialGCs = dfGCs.loc[dfGCs['ID'].isin(specialGCs)]
dfSpecialGCs = dfSpecialGCs.reset_index()

print(dfSpecialGCs)

starsX, starsY, starsZ = dfStars['xcoord'], dfStars['ycoord'], dfStars['zcoord']

fig, axs = plt.subplots(3, 3, constrained_layout=False, figsize=(10,10))

for j in [0,1,2]: #xy, xz, yz
    
    for i, row in dfSpecialGCs.iterrows(): 
    
        ax = axs[i,j]
        
        print('')
        print('')
        print(row['Name'])
        x,y,z = row['X']-8., row['Y'], row['Z']
        
        print(x,y,z)
        
        galMass = 0.
        
        #bulge, bar, disk
        for massElement in MWPotential2014:
        
        	galMass += massElement.mass(np.sqrt(x**2. + y**2.), z=z)*bovy_conversion.mass_in_msol(220., 8.)
    
        gcDist = np.sqrt(x**2. + y**2. + z**2.) 
        jacobiRadius = 3. * gcDist * (gcMasses[row['ID']] / galMass)**(1/3.)
        
        print('enclosed mass: %.02e'%(galMass))
        print('gcDist: %.03f kpc'%(gcDist))
        print('jacobiRadius: %.03f pc'%(jacobiRadius*1000.))
        
        if row['Name'] == 'omega Cen':
            ax.annotate(r'$\omega$ Cen' , xy = (0.6, 0.8), xycoords = 'axes fraction', fontsize = 16)
        else:
            ax.annotate(row['Name'] , xy = (0.6, 0.8), xycoords = 'axes fraction', fontsize = 16)
        
        if j == 0:
        
            ax.set_xlim(x-1.5, x+1.5)
            ax.set_ylim(y-1.5, y+1.5)
        
            circle1 = plt.Circle((x,y), jacobiRadius, color='r', alpha=0.2, label='Jacobi radius')
            ax.add_patch(circle1)
            
            
            ax.scatter(dfSSs['X'], dfSSs['Y'], c='C1', 
                        alpha=1.0, edgecolors='C1', marker='.', 
                        label='Forest RRL streams', linewidths=0., s=48)
            
            ax.scatter(dfOthers['X'], dfOthers['Y'], c='C2', 
                        alpha=1.0, edgecolors='C2', marker='.', 
                        label='``Other`` groupings', linewidths=0., s=48)
            
            ax.scatter(dfNewGCs['X'], dfNewGCs['Y'], c='C0', 
                        alpha=1.0, edgecolors='C0', marker='.', 
                        label='GC-like groupings', linewidths=0., s=48)
            
            ax.scatter(starsX, starsY, c='k', 
                        alpha=0.8, edgecolors='k', marker='.', 
                        label='Field RRLs', linewidths=0., s=8)
            
            ax.set_aspect('equal')
            ax.set_xticks([x-1., x-0.5, x, x+0.5, x+1.])
            ax.set_yticks([y-1., y-0.5, y, y+0.5, y+1.])
            ax.set_xlabel(r'$x_{\rm MW}$ (kpc)', fontsize=12)
            ax.set_ylabel(r'$y_{\rm MW}$ (kpc)', fontsize=12)
            
        if j == 1:
        
            ax.set_xlim(x-1.5, x+1.5)
            ax.set_ylim(z-1.5, z+1.5)
        
            circle1 = plt.Circle((x,z), jacobiRadius, color='r', alpha=0.4, label='Jacobi radius')
            ax.add_patch(circle1)
            
            ax.scatter(dfSSs['X'], dfSSs['Z'], c='C1', 
                        alpha=1.0, edgecolors='C1', marker='.', 
                        label='Forest RRL streams', linewidths=0., s=48)
            
            ax.scatter(dfOthers['X'], dfOthers['Z'], c='C2', 
                        alpha=1.0, edgecolors='C2', marker='.', 
                        label='``Other`` groupings', linewidths=0., s=48)
            
            ax.scatter(dfNewGCs['X'], dfNewGCs['Z'], c='C0', 
                        alpha=1.0, edgecolors='C0', marker='.', 
                        label='GC-like groupings', linewidths=0., s=48)
            
            ax.scatter(starsX, starsZ, c='k', 
                        alpha=0.8, edgecolors='k', marker='.', 
                        label='Field RRLs', linewidths=0., s=8)
            
            ax.set_aspect('equal')
            ax.set_xticks([x-1., x-0.5, x, x+0.5, x+1.])
            ax.set_yticks([z-1., z-0.5, z, z+0.5, z+1.])
            ax.set_xlabel(r'$x_{\rm MW}$ (kpc)', fontsize=12)
            ax.set_ylabel(r'$z_{\rm MW}$ (kpc)', fontsize=12)
        
        if j == 2:
        
            ax.set_xlim(y-1.5, y+1.5)
            ax.set_ylim(z-1.5, z+1.5)
        
            circle1 = plt.Circle((y,z), jacobiRadius, color='r', alpha=0.4, label='Jacobi radius')
            ax.add_patch(circle1)
            
            ax.scatter(dfSSs['Y'], dfSSs['Z'], c='C1', 
                        alpha=1.0, edgecolors='C1', marker='.', 
                        label='Forest RRL streams', linewidths=0., s=48)
            
            ax.scatter(dfOthers['Y'], dfOthers['Z'], c='C2', 
                        alpha=1.0, edgecolors='C2', marker='.', 
                        label='``Other`` groupings', linewidths=0., s=48)
            
            ax.scatter(dfNewGCs['Y'], dfNewGCs['Z'], c='C0', 
                        alpha=1.0, edgecolors='C0', marker='.', 
                        label='GC-like groupings', linewidths=0., s=48)
            
            ax.scatter(starsY, starsZ, c='k', 
                        alpha=0.8, edgecolors='k', marker='.', 
                        label='Field RRLs', linewidths=0., s=8)
            
            ax.set_aspect('equal')
            ax.set_xticks([y-1., y-0.5, y, y+0.5, y+1.])
            ax.set_yticks([z-1., z-0.5, z, z+0.5, z+1.])
            ax.set_xlabel(r'$y_{\rm MW}$ (kpc)', fontsize=12)
            ax.set_ylabel(r'$z_{\rm MW}$ (kpc)', fontsize=12)
          
        if i == 0 and j == 0:
            ax.legend(loc='lower right', fontsize=6)

fig.tight_layout()
plt.savefig('threeGCs.pdf', bbox_inches='tight')