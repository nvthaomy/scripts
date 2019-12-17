#========================#
# Calculate 1D histogram 
# only for one species
# kevinshen@ucsb.edu
#========================#
# Example usage:
# python 1d-histogram.py trajectory.dcd topology.pdb -axis 2 -Lx 5 -Ly 5 -Lz 10

import mdtraj, os, re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

showPlots = True
try:
  os.environ["DISPLAY"] #Detects if display is available
except KeyError:
  showPlots = False
  matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
#plt.style.use('seaborn-dark')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)

import argparse as ap
parser = ap.ArgumentParser(description="Get 1D histogram, assuming single species")

parser.add_argument('coordfile',type=str, help="trajectory file")
parser.add_argument('topfile',type=str, help="topology file")
parser.add_argument('-axis',type=int, default=0, choices=[0,1,2], help="axis to bin")
parser.add_argument('-Lx',type=float, default=10, help="Lx, default 10")
parser.add_argument('-Ly',type=float, default=10, help="Ly, default 10")
parser.add_argument('-Lz',type=float, default=10, help="Lz, default 10")
parser.add_argument('-nbins',type=int, default=100, help="Number of bins")
args = parser.parse_args()

ax = args.axis
coordfile = args.coordfile
topfile = args.topfile

print("... Loading Trajectory ...")
traj = mdtraj.load(coordfile,top=topfile)
print("... Done Loading ...")

Lx,Ly,Lz = traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2] #assuming constant box shape
box = np.array([traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2]]) #assuming constant box shape

x   = 0.5
V   = Lx*Ly*Lz
Ntot = traj.xyz.shape[1]
nbins = args.nbins

print("box: [Lx,Ly,Lz]".format(Lx, Ly, Lz))
print("V: {}".format(V))

L = box[ax]
A = np.prod(box)/L
xs = traj.xyz[:,:,ax]
xmedian = np.median(xs,1)
xs = xs - xmedian[:,None]
xs = np.ravel(xs)
dx = L/nbins

xs = np.mod(xs,L) #wrap pbc
#histA,bins = np.histogram(xs, bins=100, density=False)
#histA = histA/traj.n_frames/A/dx
hist,bins = np.histogram(xs, bins=nbins, density=True)
hist = hist*Ntot/A  #gives density per length
binmid = 0.5*(bins[1:]+bins[0:-1])


data = np.vstack([binmid,hist]).T
np.savetxt('xhist.dat',data,header='bin-midpt\tHistogram')

fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
axs.plot(binmid, hist, marker=None,ls='-',lw=0.75,mfc="None",ms=2)
#axs.legend(loc='best',prop={'size': 5})
plt.xlabel('{}'.format(['x','y','z'][ax]))
plt.ylabel('Count')
title ='{}histogram'.format(['x','y','z'][ax])
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
                
