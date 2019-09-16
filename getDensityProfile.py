#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 09:36:41 2019
@author: nvthaomy
"""
import numpy as np
import argparse, os, sys
import matplotlib.pyplot as plt
import matplotlib
"""Plotting density profile in a specific axis, currently support lammpstrj and pdb formats
   For volume fraction plot, assumes all species atoms have same size. 
	volume fraction = #atoms of interest in slab/ #all atoms in slab
"""
showPlots = True
try:
    os.environ["DISPLAY"] #Detects if display is available
except KeyError:
    showPlots = False
    matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",required=True,help="traj file, .lammpstrj or .pdb")
parser.add_argument("-x","--x",required=True, nargs = 2, type=float ,help="lower and upper bounds of box dimension in the x-direction")
parser.add_argument("-y","--y",required=True, nargs = 2, type=float ,help="lower and upper bounds of box dimension in the y-direction")
parser.add_argument("-z","--z",required=True, nargs = 2, type=float ,help="lower and upper bounds of box dimension in the z-direction")     
parser.add_argument("-ns",default = 100, type = int, help="number of slabs")
#parser.add_argument("-N", type = float, help="number of atoms in molecule of interest")
parser.add_argument("-at", nargs='+',required=True, help="atom type for molecule of interest")
parser.add_argument("-ax", required  = True, help="axis to perform slicing on, x or y or z")
parser.add_argument("-stride",default = 1, type =int, help="trajectory is read every 'stride' frames")
args = parser.parse_args()

def getindices(traj,axis):
    """get column indices for coordinate and atomtype in traj file"""
    trjFile = open(traj,'r')
    line =trjFile.readline()
    gotIndices = False
    if trajFormat == 'lammpstrj':
        while len(line):
            if "ITEM: ATOMS" in line:
                atomtype_ind = (line.split()).index('type') -2
                z_ind = (line.split()).index(axis) -2
                gotIndices = True
            if gotIndices:
                line = ""
            line =trjFile.readline()

        
    return atomtype_ind,z_ind

def plot(zs,rho,vol_frac,axis):
    z_plot = []
    nrows= 2
    ncols = 1
    plt.figure()    
    for i in range(len(zs)-1): #making z array for plotting 
        zavg = (zs[i]+zs[i+1])/2
        z_plot.append(zavg)
        
    plt.subplot(nrows,ncols,1)
    plt.plot(z_plot,rho)
    plt.ylabel('$\\rho_{monomer}$ (#/Vol)')
    plt.xlabel(axis) 
    plt.xlim(min(zs),max(zs))
    plt.title('Density Profile',loc='center')
    
    plt.subplot(nrows,ncols,2)
    plt.plot(z_plot,vol_frac)
    plt.ylabel('$\\phi_{monomer}$ (#/Ntot)')
    plt.xlabel(axis)
    plt.xlim(min(zs),max(zs))
    
    plt.subplots_adjust(hspace=0.5)
    plt.savefig('DensityProfile.png',dpi=500,transparent=True)
    plt.show()                    
    
def getDensityProfile_lammps(L,ax_ind,ns,traj,at,stride):
    "get density profile in the direction determined by ax_ind"
    axis = ['x','y','z'][ax_ind] #axis in which density profile is calculated
    L_slice = L[ax_ind] #lower and upper bounds of slicing direction 
    L1 = L[min([i for i in range(3) if i != ax_ind])] #lower and upper bounds in other two directions
    L2 = L[max([i for i in range(3) if i != ax_ind])]
    boxVol = (L_slice[1] -L_slice[0]) * (L1[1]-L1[0]) * (L2[1]-L2[0])
    zs = np.linspace(L_slice[0], L_slice[1], ns+1,endpoint = True)
    dz = (L_slice[1] -L_slice[0])/float(ns)
    slabVol = dz * (L1[1]-L1[0]) * (L2[1]-L2[0])
    atomtype_ind,z_ind = getindices(traj,axis)
    rho = np.zeros(len(zs)-1)
    N_tot = 0
    Nframe = 0
    frame = 0
    nextFrame = 1
    trjFile = open(traj,'r')
    line =trjFile.readline()
    while len(line): 
        if "ITEM: TIMESTEP" in line:
            frame  += 1
            if frame != 1 and np.sum(Natoms)!= 0:
                vol_frac_temp = rho_temp/Natoms
                vol_frac += vol_frac_temp
                rho += rho_temp
            rho_temp = np.zeros(len(zs)-1) #density list of one frame
            vol_frac_temp = np.zeros(len(zs)-1) #volume fraction list of one frame
            Natoms = np.zeros(len(zs)-1) #number of all atoms in each slab
            if frame == nextFrame:
                Nframe += 1
                readTraj = True  
                sys.stdout.write("\rReading frame {}".format(frame)) 
                sys.stdout.flush()
                nextFrame += stride
                readFrame = True
            else: 
                readFrame = False  
            readTraj = False
            #print "Frame: {}".format(frame)
        else:
            if "ITEM: ATOMS" in line:
                readTraj = True
            if readTraj and readFrame  and not "ITEM: ATOMS" in line:
                if frame == 1: #count the total number of atoms
                    N_tot += 1
                atomtype = line.split()[atomtype_ind] 
                current_z = float(line.split()[z_ind])
                if current_z < L_slice[0]: #if atom is outside of box, wrap it
                        current_z = L_slice[1]-(L_slice[0]-current_z)
                elif current_z > L_slice[1]:
                        current_z = L_slice[0]+(current_z - L_slice[1])
                bound = min(zs, key=lambda x:abs(x-current_z)) #value of coordinate in zs array that is closest to the position of atom in the slicing direction   
                if current_z <= bound:                        
                        index = np.where(zs==bound)[0][0]-1
                else:
                        index = np.where(zs==bound)[0][0]
		if atomtype in at:
                    rho_temp[index] += 1
		Natoms += 1
       
        line =trjFile.readline()
    sys.stdout.write("\nTotal number of atoms: {}".format(N_tot))
    vol_frac = vol_frac/Nframe
    rho = rho/Nframe/slabVol
    plot(zs,rho,vol_frac,axis)
    z_plot = []
    for i in range(len(zs)-1): #making z array for plotting 
        zavg = (zs[i]+zs[i+1])/2
        z_plot.append(zavg)
    data = open('densityProfile.dat','w')
    data.write('\n# r Density VolFrac')
    for i, z in enumerate(z_plot):
        data.write('\n{z} {rho} {vol_frac}'.format(z=zs[i],rho=rho[i],vol_frac=vol_frac[i]))
    return zs, rho, vol_frac   
  
def getDensityProfile_pdb(L,ax_ind,ns,traj,at,stride):
    axis = ['x','y','z'][ax_ind] 
    L_slice = L[ax_ind] #lower and upper bounds of slicing direction 
    L1 = L[min([i for i in range(3) if i != ax_ind])] #lower and upper bounds in other two directions
    L2 = L[max([i for i in range(3) if i != ax_ind])]
    boxVol = (L_slice[1] -L_slice[0]) * (L1[1]-L1[0]) * (L2[1]-L2[0])
    zs = np.linspace(L_slice[0], L_slice[1], ns+1,endpoint = True)
    dz = (L_slice[1] -L_slice[0])/float(ns)
    slabVol = dz * (L1[1]-L1[0]) * (L2[1]-L2[0])
    atomtype_ind = [12,15] #lower and upper indices of column that contains info about atom type in pdb format
    if axis == 'x':
        z_ind = [30,37] #lower and upper indices of column that contains info about coordinate in pdb format
    elif axis == 'y':
        z_ind = [38,45]
    else:
        z_ind = [46,53]
    rho = np.zeros(len(zs)-1) #density, #atoms/Vol
    N_tot = 0
    vol_frac = np.zeros(len(zs)-1) 
    readTraj = True
    Nframe = 0
    frame = 0
    nextFrame = 1
    trjFile = open(traj,'r')
    line =trjFile.readline()
    
    while len(line): 
        if "MODEL" in line:
            frame  += 1
            if frame != 1 and np.sum(Natoms)!= 0:
		rho_temp = (rho_temp)
                Natoms = (Natoms)
		vol_frac_temp = rho_temp/Natoms
                vol_frac += vol_frac_temp
		rho += rho_temp
	    rho_temp = np.zeros(len(zs)-1) #density list of one frame
	    vol_frac_temp = np.zeros(len(zs)-1) #volume fraction list of one frame
	    Natoms = np.zeros(len(zs)-1) #number of all atoms in each slab
            if frame == nextFrame:
                Nframe += 1
                readTraj = True  
                sys.stdout.write("\rReading frame {}".format(frame)) 
                sys.stdout.flush()
                nextFrame += stride
                #N_tot = np.zeros(len(zs)-1) #number of all atoms in each slice                 
            else: 
                readTraj = False
        elif ("HETATM" in line or "ATOM" in line) and readTraj:
            atomtype = line[atomtype_ind[0]:atomtype_ind[1]+1].split()[0]
            current_z = float(line[z_ind[0]:z_ind[1]+1].split()[0]) 
            if frame == 1:
                N_tot +=1             
            if current_z < L_slice[0]: #usually needed with traj from openMM. Wrap if atom is outside of box, wrap it
                current_z = L_slice[1]-(L_slice[0]-current_z)
            elif current_z > L_slice[1]:
                current_z = L_slice[0]+(current_z - L_slice[1])
            bound = min(zs, key=lambda x:abs(x-current_z)) #value of coordinate in zs array that is closest to the position of atom in the slicing direction   
            if current_z <= bound:                        
                index = np.where(zs==bound)[0][0]-1
            else:
                index = np.where(zs==bound)[0][0]
            if atomtype in at:
                rho_temp[index] += 1
            Natoms[index] +=1
          
        line =trjFile.readline()
    #vol_frac = rho/Nframe/(N_tot/boxVol*slabVol) #volume fraction, #atoms/#Ntot assuming all atoms of all species have same volume
    sys.stdout.write("\nTotal number of atoms: {}".format(N_tot))
    vol_frac = vol_frac/float(Nframe)
    rho = rho/Nframe/slabVol
    plot(zs,rho,vol_frac,axis)
    z_plot = []
    for i in range(len(zs)-1): #making z array for plotting 
        zavg = (zs[i]+zs[i+1])/2
        z_plot.append(zavg)
    data = open('densityProfile.dat','w')
    data.write('\n# r Density VolFrac')
    for i, z in enumerate(z_plot):
        data.write('\n{z} {rho} {vol_frac}'.format(z=z,rho=rho[i],vol_frac=vol_frac[i])) 
    return zs ,rho , vol_frac 
                    
#Inputs
ns = args.ns
traj = args.input
#N = args.N
stride = args.stride
trajFormat = traj.split('.')[-1]
at = args.at 
x = args.x
y = args.y
z = args.z
L = [x,y,z]
axis = ['x','y','z']
ax_ind = axis.index(args.ax) #axis = [x,y,z]
for i in L:
    if len(i) != 2:
        raise Exception('Wrong number of inputs for the lower and upper bounds of box dimensions')
    if i[0]>i[1]: #check if bounds are in the correct order
        imax = max(i)
        i[0] = min(i)
        i[1] = imax
sys.stdout.write('Getting density profile across the {axis} axis \n'.format(axis=axis[ax_ind]))
sys.stdout.write('Reading from {trajFormat} file every {stride} frame(s)\n'.format(stride = stride,trajFormat = trajFormat))
sys.stdout.write('Number of bins: {}\n'.format(ns))
if trajFormat == 'lammpstrj':
    getDensityProfile_lammps(L,ax_ind,ns,traj,at,stride)     
elif trajFormat == 'pdb':        
    getDensityProfile_pdb(L,ax_ind,ns,traj,at,stride) 
else:
    raise Exception('Only support .lammpstrj and .pdb format')
                
            
    

            
