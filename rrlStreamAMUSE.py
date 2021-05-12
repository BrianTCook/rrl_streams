#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 12:16:55 2021

@author: BrianTCook
"""

import numpy as np
import pandas as pd

#from amuse.lab import *

from astropy.coordinates import SkyCoord
import astropy.units as u
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
from galpy.potential import MWPotential2014

np.random.seed(42)
    
class streamModel:
    def __init__(self, centralMass, rrlMass, orbit, nOrbiters):

        #measured mass, phase space data for the globular cluster
        self.centralParticle = Particle()
        self.centralParticle.mass = centralMass
        self.centralParticle.x = orbit.x
        self.centralParticle.y = orbit.y
        self.centralParticle.z = orbit.z
        self.centralParticle.vx = orbit.vx
        self.centralParticle.vy = orbit.vy
        self.centralParticle.vz = orbit.vz
            
        self.enclosedGalaxyMass = np.sum( [ massElement.mass(np.sqrt(self.CentralParticle.x()**2. + self.CentralParticle.y()**2.), z=self.CentralParticle.z())*bovy_conversion.mass_in_msol(220., 8.) for massElement in MWPotential2014 ] )

        #enclosed galactic mass and jacobi radius
        self.galactocentricDist = np.sqrt(self.CentralParticle.x()**2. + self.CentralParticle.y()**2. + self.CentralParticle.z()**2.)
        self.jacobiRadius = 3. * self.galactocentricDist * (self.centralMass / self.enclosedMass)**(1/3.)
        
        self.orbitingParticles = [ Particle() for i in range(nOrbiters) ]
        
        for orbiter in self.orbitingParticles:
            
            orbiter.mass = rrlMass
            orbiter.x = orbit.x
            orbiter.y = orbit.y
            orbiter.z = orbit.z
            orbiter.vx = orbit.vx
            orbiter.vy = orbit.vy
            orbiter.vz = orbit.vz

            #need to add random theta, phi values and displace st orbiter.pos = gc.pos + \vec(jacobiradius)

            randomPhi = 2. * np.pi * np.random.rand()
            randomTheta = np.arccos(1. - 2. * np.random.rand())
            
            orbiter.x += self.jacobiRadius * np.cos(randomTheta) * np.cos(randomPhi)
            orbiter.y += self.jacobiRadius * np.cos(randomTheta) * np.sin(randomPhi)
            orbiter.z += self.jacobiRadius * np.sin(randomTheta)

def orbitGetter(dataFile, specialGCs):
    
    df = pd.read_csv(dataFile, delim_whitespace=True, names=['Name','ra','dec','dist','disterr','vlos','vloserr',
                                                        'sigma','rmax','pmra','pmdec','pmra_e','pmdec_e','pmcorr',
                                                        'plx','plx_e','not sure'])
    
    dfSub = df.loc[df['Name'].isin(specialGCs)]
    dfSub = dfSub.reset_index()
    
    sortedNames = [ '' for i in range(len(specialGCs)) ]
    
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
        sortedNames[index] = row['Name']

    return sortedNames, listOfOrbits

def solver_codes_initial_setup(galaxy_code, streamModels):
    
	'''
	will need to ask SPZ if he meant for field, orbiter to be separate in non
	Nemesis gravity solvers?
	'''

    converter_parent = nbody_system.nbody_to_si(1e11|units.MSun, 15.|units.kpc)
    converter_sub = nbody_system.nbody_to_si(0.65|units.MSun, 1.|units.parsec) 
    
    gravity = bridge.Bridge(use_threading=False)
    
    for stream in streamModels:
    
        herm = Hermite(converter_parent)
        herm.particles.add_particles(stream.centralParticle)
        #herm.particles.add_particles(stream.orbitingParticles)
        gravity.add_system(herm, (galaxy_code,))  
    
    return gravity.particles, gravity


def simulation(streamModels):

	'''
	runs three stream models
	'''

	galaxy_code = to_amuse(MWPotential2014, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
	stars_g, gravity = solver_codes_initial_setup(galaxy_code, streamModels) #stars for gravity, stars for stellar

	t_end, dt = 1000.|units.Myr, 0.01|units.Myr

	sim_times_unitless = np.arange(0., (t_end+dt).value_in(units.Myr), dt.value_in(units.Myr))
    sim_tmes = [ t|units.Myr for t in sim_times_unitless ]

	#for 3D numpy array storage
	Nsavetimes = 99
	Ntotal = len(gravity.particles)
	    
	grav_data = np.zeros((Nsavetimes, Ntotal, 7))
	energy_data = np.zeros(Nsavetimes)

	#for saving in write_set_to_file
	filename_grav = 'data_temp_grav.csv'

	attributes_grav = ('mass', 'x', 'y', 'z', 'vx', 'vy', 'vz')

	print('len(sim_times) is', len(sim_times))
	saving_flag = int(math.floor(len(sim_times)/(Nsavetimes-1)))
	print('saving_flag: %i'%(saving_flag))

	snapshot_galaxy_masses = [ 0. for i in range(Nsavetimes) ]
	snapshot_times = [ 0. for i in range(Nsavetimes) ]
	j_like_index = 0

	t0 = time.time()

	channel_from_gravity_to_framework = gravity.particles.new_channel_to(stars_g)
	channel_from_framework_to_gravity = stars_g.new_channel_to(gravity.particles)

	energy_init = gravity.particles.potential_energy() + gravity.particles.kinetic_energy()

	for j, t in enumerate(sim_times):

		if j%saving_flag == 0:

			print('j = %i'%(j))
            
			energy = gravity.particles.potential_energy() + gravity.particles.kinetic_energy()
			deltaE = energy/energy_init - 1.

			print_diagnostics(t, t0, stars_g, energy, deltaE)

			energy_data[j_like_index] = deltaE

			#gravity stuff

			io.write_set_to_file(gravity.particles, filename_grav, 'csv',
					 attribute_types = (units.MSun, units.parsec, units.parsec, units.parsec, units.kms, units.kms, units.kms),
					 attribute_names = attributes_grav)

			data_t_grav = pd.read_csv(filename_grav, names=list(attributes_grav))
			data_t_grav = data_t_grav.drop([0, 1, 2]) #removes labels units, and unit names

			data_t_grav = data_t_grav.astype(float) #strings to floats
            
			grav_data[j_like_index, :len(data_t_grav.index), :] = data_t_grav.values
			np.savetxt('PhaseSpace_%s_%s_frame_%s_LCC.ascii'%(forward_or_backward, background_str, str(j).rjust(5, '0')), data_t_grav.values)

			j_like_index += 1

		channel_from_framework_to_gravity.copy()
		gravity.evolve_model(t)
		channel_from_gravity_to_framework.copy()
    
	gravity.stop()
		
    return 1


if __name__ in '__main__':
    
    dataFile = 'gcVasiliev.txt'
    specialGCs = [ 'NGC_362', 'NGC_6626_M_28', 'Pal_5' ]
    galaxy_code = MWPotential2014
    
    sortedNames, orbits = orbitGetter(dataFile, specialGCs)
    specialGCs = sortedNames #sorted after being filled in from pandas dataframe
    
    models = [ 0. for i in range(len(specialGCs)) ]
    
    for idx, orbit in enumerate(orbits):
        
        models[idx] = streamModel(masses[idx], orbit)
        
    simulation(streamModels)
    
        