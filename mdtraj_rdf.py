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
parser.add_argument('topfile',type=str, help="topology file")
parser.add_argument('atom1', type=str, help='name of the 1st atom type')
parser.add_argument('atom2', type=str, help='name of the 2nd atom type')
parser.add_argument('-nbins',type=int, default=1000, help="Number of bins")
parser.add_argument('-rmax',type=float, default=1., help="Max distance")
parser.add_argument('-stride',type=int, default=1, help="stride")
args = parser.parse_args()
#####
coordfile = args.coordfile
topfile = args.topfile
atom1 = args.atom1
atom2 = args.atom2
nbins = args.nbins
rmax = args.rmax
stride = args.stride

print("... Loading Trajectory ...")
traj = mdtraj.load(coordfile,top=topfile,stride=stride)
top = traj.topology
print("... Done Loading ...")
Lx,Ly,Lz = traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2] #assuming constant box shape
box = np.array([traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2]]) #assuming constant box shape

V   = Lx*Ly*Lz

atoms1 = top.select("name '{}'".format(atom1))
atoms2 = top.select("name '{}'".format(atom2))

pairs = top.select_pairs(selection1=atoms1, selection2=atoms2)
r,g_r=mdtraj.compute_rdf(traj,pairs=pairs,r_range=(0,rmax),n_bins=nbins)

name='rdf_{}_{}'.format(atom1,atom2)
data = np.vstack([r, g_r]).T
np.savetxt(name+'.dat',data,header='r\tg_r')

fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
axs.plot(r, g_r, marker=None,ls='-',lw=0.75,mfc="None",ms=2)
plt.xlabel('r (nm)')
plt.ylabel('$g(r)$')
title ='rdf {} {}'.format(atom1,atom2)
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
