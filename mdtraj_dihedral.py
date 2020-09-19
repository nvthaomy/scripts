#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 21:36:10 2020

@author: nvthaomy
"""
import mdtraj,matplotlib ,os, re
import matplotlib.pyplot as plt
import argparse
import numpy as np

showPlots = True
try:
  os.environ["DISPLAY"] #Detects if display is available
except KeyError:
  showPlots = False
  matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
#plt.style.use('seaborn-dark')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)


parser = argparse.ArgumentParser()
parser.add_argument('coordfile',type=str, help="trajectory file")
parser.add_argument('topfile', help="topology file")
parser.add_argument('atoms', type=str, nargs = 4, help='name of 4 atoms involved in dihedral torsion, first atom is the center atom')
parser.add_argument('-topdat', type=str, help='if provided will make a new topology, text file of molecule pdbs (1st column) and number (2nd column) ')
parser.add_argument('-nbins',type=int, default=500, help="Number of bins")
parser.add_argument('-stride',type=int, default=1, help="stride")
parser.add_argument('-w', type=int, default=0, help='warmup frames')
parser.add_argument('-plot', type=int, nargs = "+",help='dihedral index to plot time series')

args = parser.parse_args()
#####
coordfile = args.coordfile
topfile = args.topfile
atoms = args.atoms
topdat = args.topdat
nbins = args.nbins
stride = args.stride
warmup = args.w
id4plot = args.plot

def MakeTop(topdat):
    print('Make new topology based on {}'.format(topdat))
    file = open(topdat,'r')
    lines = file.readlines()
    chainTops = []
    nMols = []
    for line in lines:
        pdb = line.split()[0]
        nMol = line.split()[1]
        nMols.append(int(nMol))
        top = mdtraj.load(pdb).topology
        chainTops.append(top)
        
    top = mdtraj.Topology()
    for ii, moltop in enumerate(chainTops):
        nummol = nMols[ii]
        for jj in range(nummol):
            atoms_in_top = []
            for c in moltop.chains:
                chain = top.add_chain()
                for r in c.residues:
                    residue = top.add_residue("custom",chain)
                    for a in enumerate(r.atoms):
                        atom = top.add_atom(a[1].name, a[1].element, residue)
                        atoms_in_top.append(atom)
                        #print(atom)
            for bond in moltop.bonds:
                bi1,bi2 = bond[0].index, bond[1].index
                top.add_bond( atoms_in_top[bi1], atoms_in_top[bi2] )
    return top

print("... Loading Trajectory ...")
traj = mdtraj.load(coordfile,top=topfile,stride=stride)
traj = traj[warmup:]
top = traj.topology
print("... Done Loading ...")
if topdat:
    top = MakeTop(topdat)
Lx,Ly,Lz = traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2] #assuming constant box shape
box = np.array([traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2]]) #assuming constant box shape
V   = Lx*Ly*Lz

firstatomIds = top.select("name '{}'".format(atoms[0]))
idsDihed = np.vstack((firstatomIds, -1 * np.ones(len(firstatomIds)), -1 * np.ones(len(firstatomIds)), -1 * np.ones(len(firstatomIds))))
idsDihed = np.array(idsDihed, dtype = int)
idsDihed = idsDihed.transpose() # list of 4 atom indices in this dihedral torsion
print('Calculating dihedral angles of {} {} {} {}'.format(*atoms))

for i, ai0 in enumerate(firstatomIds):
    ai0_current = ai0
    row = i
    for bond in top.bonds:
        if ai0_current == bond[0].index:
            if bond[1].name in atoms and not bond[1].index in idsDihed[row] :
                col = np.where(idsDihed[row] == -1)[0][0]
                idsDihed[row,col] = bond[1].index
                ai0_current = bond[1].index
        #print(np.where(idsDihed[row] == -1)[0])
        if len(np.where(idsDihed[row] == -1)[0]) == 0:
            break

# eliminate any row that does not have 4 indices by looking for -1
i = np.where(idsDihed==-1)[0]
j = np.array([k for k in range(len(firstatomIds)) if not k in i])
if len(j) > 0:
    idsDihed = idsDihed[j,:]
print('\nidsDihed \n',idsDihed)

xyz = np.zeros((idsDihed.shape[0], traj.n_frames, idsDihed.shape[1], 3))
Imprp = np.zeros((idsDihed.shape[0], traj.n_frames))
for ii, i in enumerate(idsDihed):
    xyz_tmp = traj.xyz[:,i,:]
    xyz[ii,:,:,:] = xyz_tmp

#theta = GetImprp(xyz)
theta = mdtraj.compute_dihedrals(traj,idsDihed) * 180./np.pi # convert to degrees 
theta = theta.transpose() #  nImprp x n_frames   

np.savetxt('dihed_{}_{}_{}_{}.dat'.format(*atoms), 
           np.hstack((np.arange(theta.shape[0]).reshape(theta.shape[0],1),theta)), 
           header = 'ImprpID theta(t)')

if id4plot:
    try:
        os.mkdir('dihed_{}_{}_{}_{}_plots'.format(*atoms))
    except:                                                                                   
        pass                                                                                  
    os.chdir('dihed_{}_{}_{}_{}_plots'.format(*atoms))   
    for i in id4plot:
        fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
        axs.plot(range(theta.shape[1]), theta[i], marker='o',ls='-',lw=0.75,mfc="None",ms=2)
        plt.ylabel('$\\theta$')
        plt.xlabel('frame')
        title ='dihed {}'.format(i)
        plt.title(title, loc = 'center')
        plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
    os.chdir('..')

# histogramming
theta = np.ravel(theta)
mean = np.mean(theta)
dtheta = 180./nbins
hist,bins = np.histogram(theta, bins=nbins, density=True)
binmid = 0.5*(bins[1:]+bins[0:-1])
data = np.vstack([binmid,hist]).T
np.savetxt('hist_dihed_{}_{}_{}_{}.dat'.format(*atoms),
           data,header='bin-midpt\tProbability')

fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
axs.plot(binmid, hist, marker='o',ls='-',lw=0.5,mfc="None",ms=2)
#plt.text(mean*1.3, 0.8*max(hist), '{}'.format(round(mean,3)), color='r', size=7)
#axs.vlines(mean,0,max(hist),colors = 'r',lw=1)
#axs.legend(loc='best',prop={'size': 5})
plt.xlabel('$\\theta (deg)$')
plt.ylabel('$P(\\theta)$')
#plt.xlim(0,100)
title ='hist_dihed_{}_{}_{}_{}'.format(*atoms)
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()

