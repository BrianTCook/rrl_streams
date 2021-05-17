#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 12:16:55 2021

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
    
class streamModel:
	def __init__(self, centralMass, orbit, rrlMass, nOrbiters):

		centralParticle = Particles(1)
		centralParticle.mass = centralMass|units.MSun
		centralParticle.x = orbit.x() | units.kpc
		centralParticle.y = orbit.y() | units.kpc
		centralParticle.z = orbit.z() | units.kpc
		centralParticle.vx = orbit.vx() | units.kms
		centralParticle.vy = orbit.vy() | units.kms
		centralParticle.vz = orbit.vz() | units.kms

		#measured mass, phase space data for the globular cluster
		self.centralParticle = centralParticle

		galMass = 0.|units.MSun

		#bulge, bar, disk
		for massElement in MWPotential2014:

			galMass += massElement.mass(np.sqrt(self.centralParticle.x.value_in(units.kpc)**2. + self.centralParticle.y.value_in(units.kpc)**2.), z=self.centralParticle.z.value_in(units.kpc))*bovy_conversion.mass_in_msol(220., 8.) | units.MSun

		self.enclosedGalaxyMass = galMass
		self.speed = np.sqrt(self.centralParticle.vx.value_in(units.kms)**2. + self.centralParticle.vy.value_in(units.kms)**2. + self.centralParticle.vz.value_in(units.kms)**2.) | units.kms

		#enclosed galactic mass and jacobi radius
		self.galactocentricDist = np.sqrt(self.centralParticle.x.value_in(units.kpc)**2. + self.centralParticle.y.value_in(units.kpc)**2. + self.centralParticle.z.value_in(units.kpc)**2.) | units.kpc
		self.jacobiRadius = 3. * self.galactocentricDist * (centralParticle.mass.value_in(units.MSun) / self.enclosedGalaxyMass.value_in(units.MSun))**(1/3.)
		self.crossingTime = self.galactocentricDist / self.speed

		orbitingParticles = Particles(nOrbiters)

		for orbiter in orbitingParticles:

			orbiter.mass = rrlMass | units.MSun
			orbiter.x = orbit.x() | units.kpc
			orbiter.y = orbit.y() | units.kpc
			orbiter.z = orbit.z() | units.kpc
			orbiter.vx = orbit.vx() | units.kms
			orbiter.vy = orbit.vy() | units.kms
			orbiter.vz = orbit.vz() | units.kms

			#need to add random theta, phi values and displace st orbiter.pos = gc.pos + vec(jacobiradius)

			phi = 2. * np.pi * np.random.random()
			theta = np.arccos(1. - 2. * np.random.random()) - np.pi/2.
			psi = 2. * np.pi * np.random.random()

			fCross = np.sqrt( constants.G * self.centralParticle.mass / self.jacobiRadius**(3.) )[0]
			

			thetaDot = fCross * np.cos(psi)
			phiDot = fCross * np.sin(theta) * np.sin(psi)

			orbiter.x += self.jacobiRadius * np.cos(theta) * np.cos(phi)
			orbiter.y += self.jacobiRadius * np.cos(theta) * np.sin(phi)
			orbiter.z += self.jacobiRadius * np.sin(theta)
			orbiter.vx += self.jacobiRadius * (-np.cos(theta)*np.sin(phi)*phiDot - np.sin(theta)*np.cos(phi)*thetaDot)
			orbiter.vy += self.jacobiRadius * (np.cos(theta)*np.sin(phi)*phiDot - np.sin(theta)*np.sin(phi)*thetaDot)
			orbiter.vz += self.jacobiRadius * np.cos(theta) * thetaDot

		self.orbitingParticles = orbitingParticles

def print_diagnostics(sim_time, t0, simulation_bodies, E_dyn, dE_dyn):

	print('-----------------------------------------------------')
	print('simulation time: ', sim_time)
	print('wall time: %.03f minutes'%((time.time()-t0)/60.))
	print('simulation_bodies.center_of_mass() in kpc: ', simulation_bodies.center_of_mass().value_in(units.kpc))
	print('E_dyn: ', E_dyn)
	print('dE_dyn: ', dE_dyn)
	print('-----------------------------------------------------')	

def orbitGetter(dataFile, specialGCs):
    
    df = pd.read_csv(dataFile, delim_whitespace=True, names=['Name','ra','dec','dist','disterr','vlos','vloserr',
                                                        'sigma','rmax','pmra','pmdec','pmra_e','pmdec_e','pmcorr',
                                                        'plx','plx_e','not sure'])
    
    dfSub = df.loc[df['Name'].isin(specialGCs)]
    dfSub = dfSub.reset_index()
    
    sortedNames = [ '' for i in range(len(specialGCs)) ]
    
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

	gravity = bridge.Bridge()

	for stream in streamModels:

		subCode = Huayno(converter_sub)

		subCode.particles.add_particles(stream.centralParticle)
		subCode.particles.add_particles(stream.orbitingParticles)
		gravity.add_system(subCode, (galaxy_code,))  

	return gravity.particles, gravity


def simulation(streamModels):

	'''
	runs three stream models
	'''

	galaxy_code = to_amuse(MWPotential2014, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
	stars_g, gravity = solver_codes_initial_setup(galaxy_code, streamModels) #stars for gravity, stars for stellar

	t_end, dt = 35.|units.Myr, 10000.|units.yr

	sim_times_unitless = np.arange(0., (t_end+dt).value_in(units.Myr), dt.value_in(units.Myr))
	sim_times = [ t|units.Myr for t in sim_times_unitless ]

	#for 3D numpy array storage
	Nsavetimes = 36
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
			np.savetxt('PhaseSpace_frame_%s_streams.ascii'%(str(j).rjust(5, '0')), data_t_grav.values)

			j_like_index += 1

		channel_from_framework_to_gravity.copy()
		gravity.evolve_model(t)
		channel_from_gravity_to_framework.copy()

	gravity.stop()

	return 1

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
    
    try:
        outputArr = np.loadtxt('tidalRadii.txt')
    
    except:
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
        
        np.savetxt('tidalRadii.txt', outputArr)
    
    plt.figure()
    outputArr = np.loadtxt('tidalRadii.txt')
    
    for k in range(nStreams):
        
        oldString = streamNames[k]
        newString = oldString.replace("_", " ")
        
        plt.plot(tVals, outputArr[:,k], label='%s'%(newString))
        
    plt.gca().set_yscale('log')
    
    plt.xlim(0, 35)
    plt.ylim(30, 1500)
    plt.axhline(y=1000, linestyle='--', c='k', linewidth=1)
    
    plt.gca().set_xlabel(r'$t_{\rm sim}$ [Myr]', fontsize=20)
    plt.gca().set_ylabel(r'$r_{\rm Jacobi}$ [pc]', fontsize=20)
    plt.annotate(r'$R_{\rm max}$', xy = (0.1, 0.91), xycoords = 'axes fraction', fontsize = 14)
    
    plt.gca().tick_params(labelsize='x-large')
    plt.legend(loc='upper right', fontsize=10)
    
    
    plt.savefig('jacobiRadii_appendix.pdf', bbox_inches='tight')
    
    return 1

def alignments(nStreams, nSnapshots, files, streamNames):
    
    t0 = time.time()
    outputArr = np.zeros((nSnapshots, nStreams))
    tVals = np.linspace(0., 35., 36)
    
    try:
        outputArr = np.loadtxt('medianAlignments.txt')
    
    except:
        for i, file in enumerate(files):
            
            print('')
            print('currently at snapshot %.0f Myr: '%(tVals[i]))
            print('wall time: %.02f minutes'%((time.time()-t0)/60.))
            
            df = pd.read_csv(file, sep=' ', names=['mass','x','y','z','vx','vy','vz'])
            dfStreams = np.array_split(df, nStreams)
            
            for j, dfStream in enumerate(dfStreams):
                
                vys,vzs = dfStream['vy'].tolist(), dfStream['vz'].tolist()
                
                vyBig, vzBig = vys[0], vzs[0]
                vyLittle, vzLittle = vys[1:], vzs[1:]
                
                vHat = [ vyBig, vzBig ] / np.linalg.norm( np.array([ vyBig, vzBig ]) )
                vsNormalized = [ [ vyl, vzl ] / np.linalg.norm(np.array([ vyl, vzl ])) for vyl, vzl in zip(vyLittle, vzLittle) ]
                cosThetas = [ np.dot(v, vHat) for v in vsNormalized ]
        
                outputArr[i,j] = np.median(cosThetas)
        
        np.savetxt('medianAlignments.txt', outputArr)
    
    plt.figure()

    for k in range(nStreams):
        
        oldString = streamNames[k]
        newString = oldString.replace("_", " ")
        
        plt.plot(tVals, outputArr[:,k], label='%s'%(newString))
    
    plt.ylim(0.94, 1.01)
    plt.xlim(0, 35)
    plt.axhline(y=0.98, linestyle='--', c='k', linewidth=1)
    
    plt.gca().set_xlabel(r'$t_{\rm sim}$ [Myr]', fontsize=20)
    plt.gca().set_ylabel(r'$\tilde{X}_{\theta}$', fontsize=20)
    plt.annotate(r'$\tilde{X}_{\theta, {\rm min}}$', xy = (0.1, 0.6), xycoords = 'axes fraction', fontsize = 14)
    
    plt.gca().tick_params(labelsize='x-large')
    plt.legend(loc='lower right', fontsize=10)
    
    plt.savefig('medianAlignments_appendix.pdf', bbox_inches='tight')
    
    return 1


if __name__ in '__main__':
    
    dataFile = 'gcVasiliev.txt'
    
    specialGCs = [ 'NGC_362', 'Pal_5', 'NGC_5466', 'Pal_12', 'NGC_2419' ]	
    gcMasses = {'NGC_362':3.45e5, 'Pal_5':1.39e4, 'NGC_5466':4.56e4, 'Pal_12':1.19e4, 'NGC_2419':9.81e5} #mSun
    rrlMass, nOrbiters = 0.65, 40 #mSun, intege
    
    galaxy_code = MWPotential2014
    
    sortedNames, orbits = orbitGetter(dataFile, specialGCs)
    specialGCs = sortedNames #sorted after being filled in from pandas dataframe
    
    plotting_flag = True
    simulation_flag = False

    if plotting_flag == True:
        
        nStreams, nSnapshots = 5, 36
        streamNames = sortedNames
    
        files = glob.glob('*.ascii')
        tidalRadii(nStreams, nSnapshots, files, streamNames)
        alignments(nStreams, nSnapshots, files, streamNames)

    if simulation_flag == True:

        from amuse.lab import *
        from amuse.couple import bridge
        from amuse.support import io
        
        streamModels = [ 0. for i in range(len(specialGCs)) ]
    
        for idx, (GC, orbit) in enumerate(zip(specialGCs, orbits)):
    
            stream = streamModel(gcMasses[GC], orbit, rrlMass, nOrbiters)
            streamModels[idx] = stream
            print('GC and crossing time in Myr: ', GC, stream.crossingTime.value_in(units.Myr))
    
        simulation(streamModels)
    
        
