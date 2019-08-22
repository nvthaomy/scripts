#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:29:13 2019

@author: nvthaomy
"""

from simtk import openmm, unit
from simtk.unit import *
from simtk.openmm import app
import numpy as np
import simtk.openmm.app.topology as topology

#========================================
###DEFINE ENERGY, LENGTH & MASS SCALES###
#========================================
epsilon  = 1 * kilojoules_per_mole
sigma = 0.1 * nanometer
mass = 12 * amu
N_av = 6.022140857*10**23 /mole
kb = 1.380649*10**(-23)*joules/kelvin* N_av #joules/kelvin/mol
#====================================
#1) SYSTEM DIMENSIONLESS PARAMETERS###
#====================================
#Molecules:
NP = 150
NS = 2400
N_purefluidbox = 6000
DOP = 24
NS = NS +N_purefluidbox
N = NP*DOP + NS
charge = 0.0 * elementary_charge
reduced_density = 6000/(40**3)
#Bonded potential:
l = 1 #reduced bond length
k = 2 * 3000 #reduced force constant #multiply by 2 if use force constant from lammps
#Pair potential: u = -A exp(-Br^2)
B = 1
A_PP = -100
A_SS = -100
A_PS = -130
reduced_nonbondedCutoff = 10
#wall potential
kWall = 1000 *kilocalories_per_mole /angstrom**2
buffer_fraction = 0.5
#======================
#2) Integration Options
#======================
useNVT = True
useNPT = False
wallPotentialOn = True
reduced_timestep = 0.005
reduced_temp = 1
reduced_Tdamp =100*reduced_timestep #time units
reduced_pressure = 1
reduced_Pdamp = 1000*reduced_timestep #time units
#=====================
#3) Simulation Options
#=====================
steps = 2e8
equilibrationSteps = 1e7
platform = openmm.Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}
#platform = openmm.Platform.getPlatformByName('CPU')

pdbReporter_warmup = app.pdbreporter.PDBReporter('trajectory_warmup.pdb', 10000)
dataReporter_warmup = app.statedatareporter.StateDataReporter('log_warmup.txt', 500, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, remainingTime=True,separator='\t')

pdbReporter = app.pdbreporter.PDBReporter('trajectory.pdb', 50000)
dataReporter = app.statedatareporter.StateDataReporter('log.txt', 1000, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, 
    totalEnergy=True, temperature=True, volume=True, density=True,remainingTime=True, separator='\t')
nonbondedMethod = openmm.CustomNonbondedForce.CutoffPeriodic
ewaldErrorTolerance = 0.0001
constraints = None
#rigidWater = True
#constraintTolerance = 0.000001

#=============================
###Converting to real units###
#=============================
nonbondedCutoff = reduced_nonbondedCutoff*sigma
dt = reduced_timestep* (mass*sigma*sigma/epsilon)**(1/2)
temperature = reduced_temp * epsilon/kb
friction = 1/(reduced_Tdamp) * (epsilon/(mass*sigma*sigma))**(1/2)
pressure = reduced_pressure * epsilon/(sigma**3) * N_av**-1
barostatInterval = int(reduced_Pdamp/reduced_timestep)

print ('temperature:{}'.format(temperature))
print ('pressure:{}'.format(pressure))
print ('friction:{}'.format(friction))
print ('barostat interval:{}'.format(barostatInterval))
#========================================
# Create a system and add particles to it
#========================================
system = openmm.System()
# Particles are added one at a time
# Their indices in the System will correspond with their indices in the Force objects we will add later
print ("Adding {} polymer atoms into system".format(NP*DOP))
for index in range(NP*DOP):
    system.addParticle(mass)
print ("Adding {} solvent atoms into system".format(NS))
for index in range(NS):
    system.addParticle(mass)
print("Total number of paricles in system: {}".format(system.getNumParticles()))
# Set the periodic box vectors:
number_density = reduced_density / sigma**3
volume = (N-N_purefluidbox) * (number_density ** -1)
box_edge = volume ** (1. / 3.)
buffer_l = buffer_fraction*box_edge
z_wall1 = buffer_l
z_wall2 = buffer_l + box_edge
box_z = 2*buffer_l + box_edge
box_vectors = np.diag([box_edge/angstrom for i in range(3)]) * angstroms
box_vectors[2] = [0,0,box_z/angstroms] * angstroms
system.setDefaultPeriodicBoxVectors(*box_vectors)

#==================
##CREATE TOPOLOGY##
#==================
#Topology consists of a set of Chains 
#Each Chain contains a set of Residues, 
#and each Residue contains a set of Atoms.
nmols = [NP, NS]
residues = [["Poly"], ["Sol"]]
PolymerAtomList = DOP*['P']
atomList = {"Poly":PolymerAtomList,"Sol":['S']} #list of atoms in each residue
elements = {"P":app.element.Element(200, 'Polymer', 'gP', mass),
            "S":app.element.Element(201, 'Solvent', 'gS', mass)}
def makeTop(NP,NS,DOP):
    nmols = [NP, NS]
    top = topology.Topology()
    for ispec in range(len(nmols)): #loop over each species
        for imol in range(nmols[ispec]): #loop over each molecule within species
            chain = top.addChain() #create chain for each molecule
            for ires in range( len(residues[ispec]) ): #polymer and solvent each has one residue
                resname = residues[ispec][ires]
                res = top.addResidue( resname, chain)
                atoms = atomList[resname]
                for atomInd,atomName in enumerate(atoms):
                    el = elements[atomName]
                    if atomInd > 0:
                        previousAtom = atom
                    atom = top.addAtom( atomName, el, res )
                    if atomInd > 0:
                        top.addBond(previousAtom,atom)
    return top
print ("Creating topology")
top = makeTop(NP,NS,DOP)

#=========================
#create HarmonicBondForce#
#=========================
bondedForce = openmm.HarmonicBondForce()
atomIndThusFar = 0 
for imol in range(NP): #loop over all polymer chain
    counter = 0
    for atomInd in range(atomIndThusFar,atomIndThusFar+DOP-1):
        bondedForce.addBond(atomInd,atomInd+1,l*sigma,k*epsilon/(sigma*sigma))
        counter += 1
    atomIndThusFar += counter + 1 #skip the last atom in polymer chain
system.addForce(bondedForce)

#=============================
#create custom nonbonded force
#=============================
Polymer = set()
Solvent = set()
for atom in top.atoms():
        if atom.residue.name in ['Poly']:
                Polymer.add(atom.index)
        else:
                Solvent.add(atom.index)
all_atoms = Polymer.union(Solvent)

#Polymer-polymer and Solvent-Solvent:
PP_nonbondedForce = openmm.CustomNonbondedForce('-APP*exp(-B*r^2)')
PP_nonbondedForce.setNonbondedMethod(nonbondedMethod)
PP_nonbondedForce.addGlobalParameter('B',B/(sigma*sigma)) #length^-2
PP_nonbondedForce.addGlobalParameter('APP',A_PP*epsilon) #energy/mol
PP_nonbondedForce.addInteractionGroup(Polymer,Polymer)
for i in range(system.getNumParticles()):
        PP_nonbondedForce.addParticle()
PP_nonbondedForce.setCutoffDistance(nonbondedCutoff)
system.addForce(PP_nonbondedForce)
print ("Number of particles in PP nonbonded force:{}".format(PP_nonbondedForce.getNumParticles()))
#Solvent-solvent:
SS_nonbondedForce = openmm.CustomNonbondedForce('-ASS*exp(-B*r^2)')
SS_nonbondedForce.setNonbondedMethod(nonbondedMethod)
SS_nonbondedForce.addGlobalParameter('B',B/(sigma*sigma)) #length^-2
SS_nonbondedForce.addGlobalParameter('ASS',A_SS*epsilon) #energy/mol
SS_nonbondedForce.addInteractionGroup(Solvent,Solvent)
for i in range(system.getNumParticles()):
        SS_nonbondedForce.addParticle()
SS_nonbondedForce.setCutoffDistance(nonbondedCutoff)
system.addForce(SS_nonbondedForce)
print ("Number of particles in SS nonbonded force:{}".format(SS_nonbondedForce.getNumParticles()))

#Polymer-solvent:
PS_nonbondedForce = openmm.CustomNonbondedForce('-APS*exp(-B*r^2)')
PS_nonbondedForce.setNonbondedMethod(nonbondedMethod)
PS_nonbondedForce.addGlobalParameter('B',B/(sigma*sigma)) #length^-2
PS_nonbondedForce.addGlobalParameter('APS',A_PS*epsilon) #energy/mol
PS_nonbondedForce.addInteractionGroup(Polymer,Solvent)
for i in range(system.getNumParticles()):
        PS_nonbondedForce.addParticle()
PS_nonbondedForce.setCutoffDistance(nonbondedCutoff)
system.addForce(PS_nonbondedForce)
print ("Number of particles in PS nonbonded force:{}".format(PS_nonbondedForce.getNumParticles()))

#===================
#create virtual wall
#===================
if wallPotentialOn:
	wall1 = openmm.CustomExternalForce('k/2*min(0,z-z_wall1)^2')
	wall1.addGlobalParameter('k',kWall) #energy/mol/area
	wall1.addGlobalParameter('z_wall1',z_wall1)

	wall2 = openmm.CustomExternalForce('k/2*max(0,z-z_wall2)^2')
	wall2.addGlobalParameter('k',kWall) #energy/mol/area
	wall2.addGlobalParameter('z_wall2',z_wall2)

	for i in Polymer: #only act on polymer
        	wall1.addParticle(i,[])
	for i in Polymer:
        	wall2.addParticle(i,[])
	system.addForce(wall1)
	system.addForce(wall2)

#=========================
#Set up wall force report
#=========================
class ForceReporter(object):
	def __init__(self, file, reportInterval):
		self._out = open(file, 'w')
		self._fileout = file
		self._out.close()
		self._reportInterval = reportInterval

	def __del__(self):
		pass
		#self._out.close()
	def describeNextReport(self, simulation):
		steps = self._reportInterval - simulation.currentStep%self._reportInterval
		return (steps, True, False, False, False, None) #(timesteps until next report, need positions, need velocities, need forces, need energies, wrapped positions)
	def report(self, simulation, state):
		self._out = open(self._fileout, 'a')
		positions = state.getPositions(asNumpy=True)
		f1=0.
		f2=0.
		kWall_dim = (kWall/N_av).value_in_unit(kilojoules/angstrom**2)
		
		for i in Polymer:
			z = positions[i][2].value_in_unit(angstrom)
			if z < z_wall1.value_in_unit(angstrom): 
				f1 += (kWall_dim*abs(z-z_wall1.value_in_unit(angstrom))) #force per particle = kJ/particle/angstrom
			if z > z_wall2.value_in_unit(angstrom):
				f2 += (kWall_dim*abs(z-z_wall2.value_in_unit(angstrom))) 
		#f1,f2: force exerted by wall 1 and wall 2 on solute atoms
		self._out.write('%g %g\n' %(f1, f2))
		self._out.close()
#===========================
## Prepare the Simulation ##
#===========================
if useNPT:
    system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = openmm.LangevinIntegrator(temperature, friction, dt)
simulation = app.Simulation(top,system, integrator, platform)
#simulation = app.Simulation(top,system, integrator, platform,platformProperties)
positions = (box_edge - 6*angstrom) * np.random.rand(DOP*NP,3)+ [0,0,z_wall1/angstrom+3]*angstrom
fluid_positions =[box_edge,box_edge,box_z] * np.random.rand(system.getNumParticles()-DOP*NP,3)
positions = np.append(positions,fluid_positions,axis=0)
simulation.context.setPositions(positions)
#Restart and Check point:
#to load a saved state: simulation.loadState('output.xml')
simulation.reporters.append(app.checkpointreporter.CheckpointReporter('checkpnt.chk', 5000))
#===========================
# Minimize and Equilibrate
#===========================
if useNVT:
	ensemble = "NVT"
else:
	ensemble = "NPT"

print ("Initialize {} simulation".format(ensemble))
print ("Running energy minimization")
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature*3)
#simulation.loadState('output_warmedup.xml')
simulation.reporters.append(pdbReporter_warmup)
simulation.reporters.append(dataReporter_warmup)
simulation.step(equilibrationSteps)
simulation.saveCheckpoint('checkpnt_warmedup.chk')
simulation.saveState('output_warmedup.xml')
#==========
# Simulate
#==========
print ("Running production")
simulation.reporters.append(pdbReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(ForceReporter('wallforces.txt', 100))
simulation.currentStep = 0
simulation.step(steps)
simulation.saveCheckpoint('checkpnt.chk')
simulation.saveState('output.xml')
