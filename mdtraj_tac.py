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
parser.add_argument('resrange', type=str, nargs = 2, help='resid of 1st and last residue of chain')
parser.add_argument('-fr', type=int, nargs = '+', default = [-1], help = 'frames') 
parser.add_argument('-topdat', type=str, help='if provided will make a new topology, text file of molecule pdbs (1st column) and number (2nd column) ')
parser.add_argument('-nbins',type=int, default=1000, help="Number of bins")
parser.add_argument('-stride',type=int, default=1, help="stride")
parser.add_argument('-w', type=int, default=0, help='warmup frames')

args = parser.parse_args()
#####
coordfile = args.coordfile
topfile = args.topfile
resrange= args.resrange
frames = args.fr
topdat = args.topdat
nbins = args.nbins
stride = args.stride
warmup = args.w

def MakeTop(topdat):
    file = open(topdat,'r')
    lines = file.readlines()
    chainTops = []
    nMols = []
    for line in lines:
        pdb = line.split()[0]
        nMol = line.split()[1]
        nMols.append(nMol)
        top = mdtraj.load(pdb).topology
        chainTops.append(top)
        
    top = mdtraj.Topology()
    for ii, moltop in enumerate(chainTops):
        nummol = nMols[ii]
        for jj in range(nummol):
            atoms_in_top = []
            for c in moltop.chains():
                chain = top.add_chain()
                for r in c.residues():
                    residue = top.add_residue("custom",chain)
                    for a in enumerate(r.atoms()):
                        atom = top.add_atom(a[1].name, a[1].element, residue)
                        atoms_in_top.append(atom)
                        print(atom)
            for bond in moltop.bonds():
                bi1,bi2 = bond[0].index, bond[1].index
                top.add_bond( atoms_in_top[bi1], atoms_in_top[bi2] )
    return top
def GetImprp(xyz):
    """ calculate improper torsion angle 
    xyz: nImprp x n_frames x 4 x 1 matrix of atom positions"""
    
    # get normal vectors of two planes for all impropers and all frames
    xyz1 = xyz[:,:,:3,:]
    xyz2 = xyz[:,:,1:,:]
    b = np.ones(xyz.shape[0], xyz.shape[1], 3, 1)
    
    n1 = np.linalg.solve(xyz1,b) # dimensions: nImprp x n_frames x 3 x 1
    n2 = np.linalg.solve(xyz2,b)
    n1t = np.transpose(n1,(0,1,3,2)) # dimensions: nImprp x n_frames x 1 x 3
    costheta = n1t @ n2 # dimensions: nImprp x n_frames x 1 x 1
    costheta = np.squeeze(costheta) # nImprp x n_frames
    norm = np.linalg.norm(n1,axis=(2,3)) * np.linalg.norm(n2,axis=(2,3)) # nImprp x n_frames
    costheta /= norm
    theta1 = np.arccos(costheta)
    theta2 = np.pi - theta1 # second solution
    theta = np.minimum(theta1,theta2) # take the smaller angle as solution
    theta *= 180./np.pi # convert to degrees
    return theta

def ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = v1@v2
    sinang = np.linalg.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)


def GetNewCoordinate(A,Cbb,H,Cbb_next):
    ''' get new coordinate vectors jx, jy, jz such that 
            A is the origin
            jx = Cbb-A / \Cbb-A\ + A (so that A is the origin)
            jz is normal to plane (Cbb,A,H) and is in 'opposite' direction to Cbb_next
            jy is jz x jx
            all coord vec has are unit vectors'''
    
    jx = (Cbb-A) / np.linalg.norm(Cbb-A)
    
    # find normal vector of this plane solving for x from Ax=b, which is also jz
    xyz = np.array([A,Cbb,H])
    b = np.ones(3)
    # 2 possibles direction for jz, pick one that is furthes away from Cbb_next
    n =  np.linalg.solve(xyz,b)
    jz1 = n/np.linalg.norm(n)
    jz2 = - n/np.linalg.norm(n)
    theta1 = ang(jz1, Cbb_next-A)
    theta2 = ang(jz2, Cbb_next-A)
    if np.abs(theta1) > np.abs(theta2):
        jz = jz1
    else:
        jz = jz2
    
    jy = np.cross(jz,jx)
    jy /= np.linalg.norm(jy)
    
    return jx,jy,jz

def Transform(jx,jy,jz, X):
    '''Transform X from cartesian coordinates to the new coordinates jx, jy, jz'''
    ix = np.array([1,0,0])
    iy = np.array([0,1,0])
    iz = np.array([0,0,1])   
    i_vec = [ix,iy,iz]
    j_vec = [jx,jy,jz]
    
    T = np.zeros((3,3))
    for j in range(3):
        for i in range(3):
            T[j,i] = j_vec[j] @ i_vec[i]
    return T @ X

chiral_name = 'CA'
# other 3 atoms attached to chiral atom
Cbb_name = 'CC' # C backbone attached to chiral
H_name = 'HA1' # H attached to chiral
R_name = 'CB' # atom of side group that directly attached to chiral

print("... Loading Trajectory ...")
traj = mdtraj.load(coordfile,top=topfile,stride=stride)
traj = traj[warmup:]
top = traj.topology
if topdat:
    top = MakeTop
print("... Done Loading ...")

chiral_id = top.select('name {} and resid {} to {}'.format(chiral_name,*resrange))[1:-1] # the first and last residues don't have chiral center
H_id = top.select('name {} and resid {} to {}'.format(H_name,*resrange))[1:-1]
R_id = top.select('name {} and resid {} to {}'.format(R_name,*resrange))[1:-1]
Cbb_id = top.select('name {} and resid {} to {}'.format(Cbb_name,*resrange)) # this will serve as C backbone of previous residue that bonded to a chiral atoms in chiral_id

print('... Getting tacticity from residues {} to {}'.format(resrange[0],resrange[1]))
tacticity_frames = []
fm = []
for frame in frames:
    print('---Frame {}---'.format(frame))
    orient = []
    for i,A in enumerate(chiral_id):
        Cbb = Cbb_id[i]
        Cbb_next = Cbb_id[i+1]
        H = H_id[i]
        R = R_id[i]
        
        A = traj.xyz[frame, A, :]
        Cbb = traj.xyz[frame, Cbb, :]
        Cbb_next = traj.xyz[frame, Cbb_next, :]
        H = traj.xyz[frame, H, :]
        R = traj.xyz[frame, R, :]
        #print('\nOriginal coordinates')
        #print('{} {}'.format(chiral_name, A))
        #print('{} {}'.format(Cbb_name, Cbb))
        #print('{} {}'.format(H_name, H))
        #print('{} {}'.format(R_name, R))
        
        jx,jy,jz = GetNewCoordinate(A,Cbb,H, Cbb_next) # the new coordinates will be that A, Cbb, H are on the xy plane (z component = 0)
        
        # after this, A, Cbb, and H will have the same z coordinate, just remove it for convenience
        At = Transform(jx,jy,jz, A) 
        Cbbt = Transform(jx,jy,jz, Cbb)
        Ht = Transform(jx,jy,jz, H)
        Rt = Transform(jx,jy,jz, R)
        z = At[2]
        At[2] -= z
        Cbbt[2] -= z
        Ht[2] -= z
        Rt[2] -= z
        
        #print('\nNew coordinates {} {} {}'.format(jx,jy,jz))
        #print('{} {}'.format(chiral_name, At))
        #print('{} {}'.format(Cbb_name, Cbbt))
        #print('{} {}'.format(H_name, Ht))
        #print('{} {}'.format(R_name, Rt))
        

        # take projection of Rt on the new xy plane (substract the z component)
        Rtp = Rt - [0,0,Rt[2]]
        fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,3])
        plt.plot(At[0],At[1],marker = 'o', label=chiral_name)
        plt.plot(Cbbt[0],Cbbt[1],marker = 'o', label=Cbb_name)
        plt.plot(Ht[0],Ht[1],marker = 'o', label=H_name)
        plt.plot(Rt[0],Rt[1],marker = 'o', label=R_name)
        plt.legend(loc='best',prop={'size':5})
        
        x = [At[0],At[1],Cbbt[0],Cbbt[1],Ht[0],Ht[1],Rt[0],Rt[1]]
#        plt.show()
        # get orientation of 3 points Ht, Cbbt, and Rtp
        X = np.cross(Ht-Cbbt,Rtp-Cbbt)[-1]

        if X > 0.0:
            orient.append(1) # clockwise
        elif X < 0.0:
            orient.append(-1) # counterclockwise
        else:
            orient.append(0) # 3 points are on the same line 

    orient = np.array(orient)
    print('orient {}'.format(orient))
    #get tacticity from orientation
    tac_tmp = orient[:-1] * orient[1:]
    tac = []
    for a in tac_tmp:
        if a > 0.:
            tac.append('m')
        elif a < 0.:
            tac.append('r')   
        else:
            tac.append('nan')
    fm.append(tac.count('m')/len(tac))    
    tacticity_frames.append(tac)           
    print('Tacticity: {}'.format(tac))
    print('Meso fraction {}\n'.format(fm))
