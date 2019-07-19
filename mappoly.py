#!/usr/bin/env python

import sim
import os, sys, pickleTraj
import cgmodel as cg
import numpy as np
import argparse as ap
#import mdtraj as md
"""map AA traj to CG traj
    pdb file must be generated by ambpdb -p <parm> -c <crd> > <pdb>"""
def readPdb(pdb):
    """get number of atoms per residue (monormer)
        number of Na+ and PAA"""
    AAatomsPerMon = []
    with open(pdb,'r') as infile:
        s = infile.readline()
        MonList = []
        numNa = 0
        atomNum = 0
        resNum_old = 1 #if first monomer in pdb belongs to PAA
        res_old = False
        while len(s):
            if s.startswith('ATOM'):
                resNum = int(s.split()[4])
                res = s.split()[3]
                if not resNum == resNum_old:
                    AAatomsPerMon.append(atomNum)
                    atomNum = 0
                    if res_old:
                        MonList.append(res_old)
                atomNum += 1 
                resNum_old = resNum
                res_old = res
            elif s.startswith('TER'):
                AAatomsPerMon.append(atomNum)
                MonList.append(res_old)
                break
            s = infile.readline() 
    with open(pdb,'r') as infile:
        s = infile.readline()
        numNa = 0
        numPAA = 0
        PAA = False
        while len(s):
            if s.startswith('ATOM'):
                if s.split()[3].startswith('Na'):
                    numNa += 1
                    PAA = False
                elif s.split()[3].startswith('A'):
                    PAA = True
                else:
                    PAA = False
            elif s.split()[0]=='TER' and PAA:
                numPAA += 1
            s = infile.readline()
    return AAatomsPerMon,MonList,numNa,numPAA

parser = ap.ArgumentParser(description="mapper for PAA")
parser.add_argument('file', type=str, help = "trajectory.nc or .lammpstrj filename")
#parser.add_argument('--nPAA', action='store', required = True, 
                  #  type=int, help = "Number of PAA")
parser.add_argument('--nSalt', action='store', type=int, 
                    default=0, help = "Number of additional salt pairs, default 0")
#parser.add_argument('--f',action='store', required = True, type=float,
                    #help = "Degree of ionization")
#parser.add_argument('--N', type = int, required = True,
                    #help = 'Degree of polymerization')
parser.add_argument('-pdb', required = True, help = 'pdb file')
parser.add_argument('-p', required = True, help = 'topology file')
args = parser.parse_args()


print("\n ===== Setting Up System =====")
#N = args.N
#f = args.f
pdb = args.pdb
nCl = args.nSalt
AAatomsPerMon,MonList,nNa,nPAA = readPdb(pdb)
nNa += args.nSalt
AAatomsPerMon = np.array(AAatomsPerMon) # which atoms in PAA map to each CG monomer, 9 atoms for head and tail bead 
print 'AAatomsPerMon'
print AAatomsPerMon
print '# Na'
print nNa
print "# PAA"
print nPAA
print 'Mon List'
print MonList

NCGMon = len(AAatomsPerMon)*nPAA # number of coarse-monomers
NAAMon = AAatomsPerMon.sum() * nPAA # number of all-atom monomers

#assign atom id to each CG bead
atomTypePAA = []
for i,res in enumerate(MonList):
    if res == 'AHP' or res == 'AP' or res == 'ATP':
        atomTypePAA.append(0)
    else:
        atomTypePAA.append(1)
atomTypeNa = [2]
atomTypeCl = [3]
print 'atomTypePAA'
print atomTypePAA


InTraj = os.path.abspath(sys.argv[1])
if '.nc' in InTraj:
    t = md.load(InTraj,top = args.p)
    lammpsTraj = InTraj.split('.nc')[0]+'.lammpstrj'
    t.save(lammpsTraj)
    InTraj = os.path.abspath(lammpsTraj)
    
# get output file name
OutTraj = InTraj.split('.lammpstrj')[0] + '_mapped.lammpstrj.gz'

# create system object
#Sys = cg.createsys(NMol=[nPAA,nNa])


# ===== create mapped object =====
print("\n ===== Creating Index Mapper =====")
# index mapping for PAA
Map = sim.atommap.PosMap()
CGmonThusFar = 0
AAatomsThusFar = 0
for ii in range(nPAA): #loop over PAA molecules
    for jj in range(len(AAatomsPerMon)): #loop over CG monomers in PAA
        Atoms1 = range(AAatomsThusFar, AAatomsThusFar+AAatomsPerMon[jj])
        Atom2 = CGmonThusFar
        this_Map = sim.atommap.AtomMap(Atoms1 = Atoms1, Atom2 = Atom2) #Atoms1=AA atoms to be mapped to Atom2, Atom2 = CG bead
        Map += [this_Map]
        CGmonThusFar += 1
        AAatomsThusFar += AAatomsPerMon[jj]

# index mapping for Na
for jj in range(nNa):
    Atoms1 = range(AAatomsThusFar, AAatomsThusFar+1)
    Atom2  = CGmonThusFar
    this_Map = sim.atommap.AtomMap(Atoms1 = Atoms1, Atom2 = Atom2)
    Map += [this_Map]
    CGmonThusFar += 1
    AAatomsThusFar += 1

# index mapping for Cl
for jj in range(nCl):
    Atoms1 = range(AAatomsThusFar, AAatomsThusFar+1)
    Atom2  = CGmonThusFar
    this_Map = sim.atommap.AtomMap(Atoms1 = Atoms1, Atom2 = Atom2)
    Map += [this_Map]
    CGmonThusFar += 1
    AAatomsThusFar += 1


# ===== create atomtypes for CG system
AtomTypes = atomTypePAA * nPAA + atomTypeNa * nNa + atomTypeCl * nCl

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

