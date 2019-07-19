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
NP = 8
NS = 1808
DOP = 24
N = NP*DOP + NS
charge = 0.0 * elementary_charge
reduced_density = 2000/27.7**3
#Bonded potential:
l = 1 #reduced bond length
k = 2 * 3000 #reduced force constant #multiply by 2 if use force constant from lammps
#Pair potential: u = -A exp(-Br^2)
B = 1
A_PP = -100
A_SS = -100
A_PS = -150
reduced_nonbondedCutoff = 10
#======================
#2) Integration Options
#======================
useNVT = True
useNPT = False
reduced_timestep = 0.005
reduced_temp = 1
reduced_Tdamp = 0.1 #time units
reduced_pressure = 1
reduced_Pdamp = 0.1 #time units
#=====================
#3) Simulation Options
#=====================
steps = 100000000
equilibrationSteps = 0 
platform = openmm.Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}
#platform = openmm.Platform.getPlatformByName('CPU')
pdbReporter = app.pdbreporter.PDBReporter('trajectory2.pdb', 10000)
dataReporter = app.statedatareporter.StateDataReporter('log2.txt', 500, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, 
    totalEnergy=True, temperature=True, volume=True, density=True, remainingTime=True,separator='\t')
dataReporter_warmup= app.statedatareporter.StateDataReporter('log_warmup.txt', 100, totalSteps=equilibrationSteps,
    step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, remainingTime=True,separator='\t')

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
volume = N * (number_density ** -1)
box_edge = volume ** (1. / 3.)
box_vectors = np.diag([box_edge/angstrom for i in range(3)]) * angstroms
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
#gaussianFunc = '-A*exp(-B*r^2)'
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

#force.setUseSwitchingFunction(True) # use a smooth switching function to avoid force discontinuities at cutoff
#force.setSwitchingDistance(0.8*nonbondedCutoff)

#===========================
## Prepare the Simulation ##
#===========================
if useNPT:
    system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = openmm.LangevinIntegrator(temperature, friction, dt)
#simulation = app.Simulation(top,system, integrator, platform)
simulation = app.Simulation(top,system, integrator, platform,platformProperties)
positions = box_edge * np.random.rand(N,3)
simulation.context.setPositions(positions)
#Restart and Check point:
#to load a saved state: simulation.loadState('output.xml')
simulation.reporters.append(app.checkpointreporter.CheckpointReporter('checkpnt.chk', 1000))
#===========================
# Minimize and Equilibrate
#===========================
if useNVT:
	ensemble = "NVT"
else:
	ensemble = "NPT"

print ("Initialize {} simulation".format(ensemble))
print ("Running energy minimization")
#simulation.minimizeEnergy(maxIterations=1000)
#simulation.context.setVelocitiesToTemperature(temperature*3)
#print ("Running equilibration")
#simulation.reporters.append(dataReporter_warmup)
#simulation.step(equilibrationSteps)
#simulation.saveState('output_warmup.xml')
simulation.loadState('output.xml')
#==========
# Simulate
#==========
print ("Running production")
simulation.reporters.append(pdbReporter)
simulation.reporters.append(dataReporter)
simulation.currentStep = 0
simulation.step(steps)
simulation.saveState('output.xml')
