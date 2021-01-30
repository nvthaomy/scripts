# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:29:13 2019

@author: My Nguyen 

Simulating Gaussian polymer solution in NVT/NPT ensemble
With or without the external potential

Flow of the script

Provide inputs in:
    0) define energy, length, mass scales
    1) dimensionless parameters defining system topology and interactions
       d) only if apply an sinusoidal external potential on the polymers, set Uext to 0 to disable this
    2) integration options
    3) set up reports

4)Convert dimensionless parameters to real units based on the scales provided in 0) since openmm only handle real units
5)-9) Create system and interactions
10)-12) Simulation
"""
#import mdtraj as md
from simtk import openmm, unit
from simtk.unit import *
from simtk.openmm import app
import numpy as np
import simtk.openmm.app.topology as topology

#========================================
#0) DEFINE ENERGY, LENGTH & MASS SCALES###
#========================================
epsilon  = 1. * kilojoules_per_mole
sigma = 1. * nanometer
mass = 12 * amu
N_av = 6.022140857*10**23 /mole
kb = 1.380649*10**(-23)*joules/kelvin* N_av #joules/kelvin/mol
#====================================
#1) SYSTEM DIMENSIONLESS PARAMETERS###
#====================================
#a)Molecules:
#number of polymers and degree of polymerization
NP = 20
DOP = 24
#number of solvents
NS = 600
N = NP*DOP + NS

charge = 0.0 * elementary_charge
#density in (length)**-3
reduced_density = 6000/(40**3)

#b)Bonded potential: k/2(r-l)**2
#reduced equilibrium bond length
l = 1 
#reduced force constant, note: k = 2*k_lammps 
k = 2 * 3000 

#c)Pair potential: u = -A exp(-Br^2)
B = 1
#interaction energy scale for each pair type
#polymer-polymer
A_PP = -100 
#solvent-solvent
A_SS = -100
#polymer-solvent
A_PS = -130
reduced_nonbondedCutoff = 10
nonbondedMethod = openmm.CustomNonbondedForce.CutoffPeriodic
#whether to add correction to the energy beyond nonbonded cutoff
UseLongRangeCorrection = False

#d) External potential: applying sinusoidal external potential on polymer centroid, set Uext = 0 to disable 
mapping = 6 #apply the external potential on the centroid of this many polymer monomers
Uext = 0 * mapping #amplitude of sinusoidal potential applied on the centroids
Nperiod = 1
axis  = 0
reduced_planeLoc = 0 
reduced_L = 80. #box length along axis where potential is applied

#======================
#2) Integration Options
#======================
useNVT = True
useNPT = False
reduced_timestep = 0.005
reduced_temp = 1
reduced_Tdamp = 100*reduced_timestep #time units
reduced_pressure = 1
reduced_Pdamp = 1000*reduced_timestep #time units

steps = 2e7
equilibrationSteps = 1e7
#if platform is not set, openmm will try to select the fastest available Platform
platform = None
platformProperties = None
#platform = openmm.Platform.getPlatformByName('CPU')
#platformProperties = {'Precision': 'mixed'}
#platform = openmm.Platform.getPlatformByName('CUDA')

#=====================
#3) Reports
#=====================
#set up data and trajectory reports:
traj = 'trajectory.dcd'
#if use mdtraj's reporter:
#dcdReporter = mdtraj.reporters.DCDReporter(traj, 5000)

#if don't have mdtraj, can use openmm's reporter to write out trajectory in pdb or dcd format
#dcdReporter = app.dcdreporter.DCDReporter(traj,10000,enforcePeriodicBox=True)
pdbReporter = app.pdbreporter.PDBReporter( traj, 5000)

dataReporter = app.statedatareporter.StateDataReporter('log.txt', 1000, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, 
    totalEnergy=True, temperature=True, volume=True, density=True, remainingTime=True,separator='\t')

dcdReporter_warmup = app.dcdreporter.DCDReporter('trajectory_warmup.dcd',10000,enforcePeriodicBox=True)
dataReporter_warmup= app.statedatareporter.StateDataReporter('log_warmup.txt', 1000, totalSteps=equilibrationSteps,
    step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, remainingTime=True,separator='\t')

#=============================
#4) Converting to real units###
#=============================
planeLoc =  reduced_planeLoc*sigma
nonbondedCutoff = reduced_nonbondedCutoff*sigma
L = reduced_L*sigma
dt = reduced_timestep* (mass*sigma*sigma/epsilon)**(1/2)
temperature = reduced_temp * epsilon/kb
friction = 1/(reduced_Tdamp) * (epsilon/(mass*sigma*sigma))**(1/2)
pressure = reduced_pressure * epsilon/(sigma**3) * N_av**-1
barostatInterval = int(reduced_Pdamp/reduced_timestep)
print ('\nParameters in real units')
print ('temperature:{}'.format(temperature))
print ('pressure:{}'.format(pressure))
print ('friction:{}'.format(friction))
print ('barostat interval:{}'.format(barostatInterval))
print ('nonbondedCutoff: {}'.format(nonbondedCutoff))

#========================================
#5) Create a system and add particles to it
#========================================
system = openmm.System()
# Particles are added one at a time
# Their indices in the System will correspond with their indices in the Force objects we will add later
print ("\nAdding {} polymer atoms into system".format(NP*DOP))
for index in range(NP*DOP):
    system.addParticle(mass)
print ("Adding {} solvent atoms into system".format(NS))
for index in range(NS):
    system.addParticle(mass)
print("Total number of paricles in system: {}".format(system.getNumParticles()))
# Set the periodic box vectors:
number_density = reduced_density / sigma**3
volume = N * (number_density ** -1)
if Uext == 0: #use cubic box if no external potential is applied
	box_edge = [volume ** (1. / 3.)] * 3
else: #use rectangular box where L is the length of box in the direction of the external potential
	box_edge_short = (volume/L)**(1./2.)
	box_edge = [box_edge_short]*3
	box_edge[axis] = L
print ('Box dimensions {}'.format(box_edge))
box_vectors = np.diag([edge/angstrom for edge in box_edge]) * angstroms
system.setDefaultPeriodicBoxVectors(*box_vectors)

#==================
#6) Create topology
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
print ("\nCreating topology")
top = makeTop(NP,NS,DOP)

#===========================
#7) create HarmonicBondForce
#===========================
bondedForce = openmm.HarmonicBondForce()
atomIndThusFar = 0 
for imol in range(NP): #loop over all polymer chain
    counter = 0
    for atomInd in range(atomIndThusFar,atomIndThusFar+DOP-1):
        bondedForce.addBond(atomInd,atomInd+1,l*sigma,k*epsilon/(sigma*sigma))
        counter += 1
    atomIndThusFar += counter + 1 #skip the last atom in polymer chain
system.addForce(bondedForce)

#================================
#8) create custom nonbonded force
#================================
#gaussianFunc = '-A*exp(-B*r^2)'
#for each pair interaction type, need to add ALL atoms in the system to the force object, only paricles in the InteractionGroup will interact with each other
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
PP_nonbondedForce.setUseLongRangeCorrection(UseLongRangeCorrection)
system.addForce(PP_nonbondedForce)
print ("\nNumber of particles in PP nonbonded force:{}".format(PP_nonbondedForce.getNumParticles()))
#Solvent-solvent:
SS_nonbondedForce = openmm.CustomNonbondedForce('-ASS*exp(-B*r^2)')
SS_nonbondedForce.setNonbondedMethod(nonbondedMethod)
SS_nonbondedForce.addGlobalParameter('B',B/(sigma*sigma)) #length^-2
SS_nonbondedForce.addGlobalParameter('ASS',A_SS*epsilon) #energy/mol
SS_nonbondedForce.addInteractionGroup(Solvent,Solvent)
for i in range(system.getNumParticles()):
        SS_nonbondedForce.addParticle()
SS_nonbondedForce.setCutoffDistance(nonbondedCutoff)
SS_nonbondedForce.setUseLongRangeCorrection(UseLongRangeCorrection)
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
PS_nonbondedForce.setUseLongRangeCorrection(UseLongRangeCorrection)
system.addForce(PS_nonbondedForce)
print ("Number of particles in PS nonbonded force:{}".format(PS_nonbondedForce.getNumParticles()))

#force.setUseSwitchingFunction(True) # use a smooth switching function to avoid force discontinuities at cutoff
#force.setSwitchingDistance(0.8*nonbondedCutoff)

#=============================
#9) create external potential
#=============================
external={"U":Uext*epsilon,"NPeriod":Nperiod,"axis":axis ,"planeLoc":planeLoc}
direction=['x','y','z']
ax = external["axis"]
#atomsInExtField = [elementMap[atomname]]
if external["U"] > 0.0 * kilojoules_per_mole:
        print('Creating sinusoidal external potential in the {} direction'.format(direction[axis]))
        energy_function = 'U*sin(2*pi*NPeriod*({axis}1-r0)/L)'.format(axis=direction[ax])
        fExt = openmm.CustomCentroidBondForce(1,energy_function)
        fExt.addGlobalParameter("U", external["U"])
        fExt.addGlobalParameter("NPeriod", external["NPeriod"])
        fExt.addGlobalParameter("pi",np.pi)
        fExt.addGlobalParameter("r0",external["planeLoc"])
        fExt.addGlobalParameter("L",box_edge[ax])
        atomThusFar = 0
        for i in range(int(NP*DOP/mapping)): #looping through CG beads
            fExt.addGroup(range(atomThusFar,atomThusFar+mapping)) #assuming the first NP*DOP atoms are polymer atoms and atom index increases along chain
            fExt.addBond([i], [])
            atomThusFar += mapping
        system.addForce(fExt)

#==============================
#10) Prepare the Simulation ##
#==============================
if useNPT:
    system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = openmm.LangevinIntegrator(temperature, friction, dt)
#use platform if specified, otherwise let omm decide
if platform:
    if platformProperties:
        simulation = app.Simulation(top,system, integrator, platform,platformProperties)
    else:
        simulation = app.Simulation(top,system, integrator, platform)
else:
    simulation = app.Simulation(top,system, integrator)

#initialize positions
positions = [box_edge[0]/angstrom,box_edge[1]/angstrom,box_edge[2]/angstrom]*angstrom * np.random.rand(N,3)
simulation.context.setPositions(positions)
#write initial positions to pdb
initialpdb = 'trajectory_init.pdb'
app.PDBFile.writeModel(simulation.topology, positions, open(initialpdb,'w'))

#Restart and Check point:
#to load a saved state: simulation.loadState('output.xml') or simulation.loadCheckpoint('checkpnt.chk')
simulation.reporters.append(app.checkpointreporter.CheckpointReporter('checkpnt.chk', 1000))

#============================
#11) Minimize and Equilibrate
#============================
if useNVT:
	ensemble = "NVT"
else:
	ensemble = "NPT"

print ("\nInitialize {} simulation".format(ensemble))
print ("Running energy minimization")
simulation.minimizeEnergy(maxIterations=1000)
simulation.context.setVelocitiesToTemperature(temperature*3)
print ("Running equilibration")
simulation.reporters.append(dataReporter_warmup)
simulation.reporters.append(dcdReporter_warmup)
simulation.step(equilibrationSteps)
simulation.saveState('output_warmup.xml')

#=============
#12) Simulate
#=============
print ("Running production")
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.currentStep = 0
simulation.step(steps)
simulation.saveState('output.xml')

#t = md.load(traj,top=initialpdb)
#trajOut = traj.split('.')[0] + '.pdb'
#t.save(trajOut)
#trajOut = traj.split('.')[0] + '.lammpstrj'
#t.save(trajOut)
