#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 12:16:55 2021

@author: BrianTCook
"""

import numpy as np
import pandas as pd

from astropy.coordinates import SkyCoord
import astropy.units as u
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
    
class streamModel:
    def __init__(self, centralMass, orbit):

        #measured mass, phase space data for the globular cluster
        self.centralMass = centralMass
        self.x = orbit.x
        self.y = orbit.y
        self.z = orbit.z
        self.vx = orbit.vx
        self.vy = orbit.vy
        self.vz = orbit.vz
        
        self.enclosedMass = np.sum( [ massElement.mass(np.sqrt(self.x()**2. + self.y()**2.), z=self.z()) for massElement in MWPotential2014 ] )

        #enclosed galactic mass and jacobi radius
        self.galactocentricDist = np.sqrt(self.x()**2. + self.y()**2. + self.z()**2.)
        self.jacobiRadius = 3. * self.galactocentricDist * (self.centralMass / self.enclosedMass)**(1/3.)

def orbitGetter(dataFile, specialGCs):
    
    df = pd.read_csv(dataFile, delim_whitespace=True, names=['Name','ra','dec','dist','disterr','vlos','vloserr',
                                                        'sigma','rmax','pmra','pmdec','pmra_e','pmdec_e','pmcorr',
                                                        'plx','plx_e','not sure'])
    
    dfSub = df.loc[df['Name'].isin(specialGCs)]
    dfSub = dfSub.reset_index()
    
    print('')
    print('')
    print(dfSub)
    
    listOfOrbits = [ '' for i in range(len(dfSub.index)) ]
    
    for index, row in dfSub.iterrows():
        
        c = SkyCoord(ra=row['ra']*u.deg,dec=row['dec']*u.deg,distance=row['dist']*u.kpc,
                pm_ra_cosdec=-1.*row['pmra']*np.cos(row['dec'])*u.mas/u.yr,
                pm_dec=row['pmdec']*u.mas/u.yr,
                radial_velocity=row['vlos']*u.km/u.s)

        listOfOrbits[index] = Orbit(c, ro=8.,vo=220.,solarmotion='hogg')

    return listOfOrbits

if __name__ in '__main__':
    
    dataFile = 'gcVasiliev.txt'
    specialGCs = [ 'NGC_362', 'NGC_5466', 'NGC_6626_M_28', 'Pal_5' ]
    
    orbits = orbitGetter(dataFile, specialGCs)
    
    print(len(orbits))
    
    for idx, (GC, orbit) in enumerate(zip(specialGCs, orbits)):
        
        gcModel = streamModel(1e5, orbit)
        print(GC)
        print(gcModel.enclosedMass)
    
    