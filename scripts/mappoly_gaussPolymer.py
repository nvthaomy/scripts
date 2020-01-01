#!/usr/bin/env python

import sim
import os, sys, pickleTraj
import cgmodel as cg
import numpy as np
import argparse as ap
#import mdtraj as md
"""map AA traj to CG traj
    pdb file must be generated by ambpdb -p <parm> -c <crd> > <pdb>"""

parser = ap.ArgumentParser(description="mapper for Gaussian polymer")
parser.add_argument('file', type=str, help = ".lammpstrj filename")
parser.add_argument('-N', required = True, type = int, help = 'DOP')
parser.add_argument('-np', required = True, type = int, help = 'number of polymer chain')
parser.add_argument('-m', default = 1, type = int, help = 'mapping ratio')
args = parser.parse_args()


print("\n ===== Setting Up System =====")
N = args.N
np = args.np
AAatomsPerMon = [1]*N

print 'AAatomsPerMon'
print AAatomsPerMon
print "# polymer"
print np

NCGMon = len(AAatomsPerMon)/args.m*np # number of coarse-monomers
NAAMon = sum(AAatomsPerMon) * np # number of all-atom monomers
print 'Number of CG monomers'
print NCGMon
#assign atom id to each CG bead
atomType = [1]+(N/args.m-2)*[2]+[1]
print 'atomType'
print atomType

InTraj = os.path.abspath(sys.argv[1])
# get output file name
OutTraj = InTraj.split('.lammpstrj')[0] + '_mapped.lammpstrj.gz'

# ===== create mapped object =====
print("\n ===== Creating Index Mapper =====")
# index mapping for PAA
Map = sim.atommap.PosMap()
CGmonThusFar = 0
AAatomsThusFar = 0
atomsPerCGbead = 0
for ii in range(np): #loop over polymer chain
    ind = 0
    for jj in range(len(AAatomsPerMon)/args.m): #loop over CG monomers in 1 chain
		for kk in range(args.m):
			atomsPerCGbead += AAatomsPerMon[ind]
			ind += 1
        	Atoms1 = range(AAatomsThusFar, AAatomsThusFar+ atomsPerCGbead)
        	Atom2 = CGmonThusFar
        	this_Map = sim.atommap.AtomMap(Atoms1 = Atoms1, Atom2 = Atom2) #Atoms1=AA atoms to be mapped to Atom2, Atom2 = CG bead
        	Map += [this_Map]
        	CGmonThusFar += 1
        	AAatomsThusFar += atomsPerCGbead
		atomsPerCGbead = 0
# ===== create atomtypes for CG system
AtomTypes = atomType * np

# ===== read AA traj =====
print("\n ===== Reading AA Traj =====")
Trj = pickleTraj(InTraj)
BoxL = Trj.FrameData['BoxL']

# ===== write out new mapped traj =====
print("\n ===== Mapping and Writing Trajectory =====")
MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes, BoxL = BoxL)

# ===== convert to lammps =====
print("\n ===== Converting to LAMMPS =====")
sim.traj.base.Convert(MappedTrj, sim.traj.LammpsWrite, FileName = OutTraj, Verbose = True)
