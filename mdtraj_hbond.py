from scipy import stats
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
parser.add_argument('-a1', type=str, nargs="+", help='names of oxygens in the 1st group')
parser.add_argument('-a2', type=str, nargs="+", help='names of oxygens in the 2nd group')
parser.add_argument('-rmax',type=float, default=0.35, help="Max distance in nm between two oxygens to be considered as H-bond")
parser.add_argument('-amin', type=float, default=130., help="Min angle in degree of O(acceptor)-H(donor)-O(donor) to be considered as H-bond")
parser.add_argument('-stride',type=int, default=1, help="stride")
parser.add_argument('-w', type=int, default=0, help='warmup frames')
parser.add_argument('-ext', type=str, required = True,help='extension to output files')
args = parser.parse_args()
#####
coordfile = args.coordfile
topfile = args.topfile
atom1 = args.a1
atom2 = args.a2
ext = args.ext
rmax = args.rmax
amin = args.amin
stride = args.stride
warmup = args.w
print("... Loading Trajectory ...")
traj = mdtraj.load(coordfile,top=topfile,stride=stride)
traj = traj[warmup:]
top = traj.topology
print("... Done Loading ...")
Lx,Ly,Lz = traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2] #assuming constant box shape
box = np.array([traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2]]) #assuming constant box shape
V   = Lx*Ly*Lz

O1 = []
O2 = []
for i in atom1:
    O1.extend(top.select("name '{}'".format(i)))
for i in atom2:
    O2.extend(top.select("name '{}'".format(i)))
pairs = top.select_pairs(selection1=O1, selection2=O2)
print('{} total OO pairs'.format(len(pairs)))

# distance criterion
print('Calculate O O distance')
r = mdtraj.compute_distances(traj,pairs) # frames x pairs
distanceBool = np.zeros_like(r, dtype=bool) # all False
print('Evaluate distance criterion, OO distance must be within {} nm'.format(rmax))
i = np.where(r <= rmax)
distanceBool[i] = True
# eliminate pairs that never in contact
distanceBool_tmp = np.sum(distanceBool,axis=0) # size = pairs
i = np.where(distanceBool_tmp == True)[0]
distanceBool = distanceBool[:,i] # frames x pairs
pairs = pairs[i]
print('{} OO pairs satisfy distance criterion'.format(len(pairs)))
# trim indices of O1 and O2 to include only atoms in the current pairs array
O1 = np.intersect1d(O1,pairs.flatten())
O2 = np.intersect1d(O2,pairs.flatten()) 

triplets = []
OOpairs = []
# OHO angle criterion
# case 1: O1 is O donor: O1-H ... O2
print('Get possible triplets if oxygens in group 1 are hbond donors')
O1H = []
OH_tmp1 = [[bond[0].index,bond[1].index] for bond in top.bonds if bond[0].index in O1 and bond[1].element.name=='hydrogen']
OH_tmp2 = [[bond[1].index,bond[0].index] for bond in top.bonds if bond[1].index in O1 and bond[0].element.name=='hydrogen'] 
O1H.extend( OH_tmp1)
O1H.extend(OH_tmp2)
for i in range(len(O1H)):
   for j in range(len(O2)):
       tmp = O1H[i].copy()
       tmp.extend([O2[j]])
       tmp = np.array(tmp,dtype=int)
       triplets.append(tmp)
       OOpairs.append([O1H[i][0],O2[j]])
# case 2: O2 is O donor: O2-H ... O1
print('Get possible triplets if oxygens in group 2 are hbond donors')
O2H = []
OH_tmp1 = [[bond[0].index,bond[1].index] for bond in top.bonds if bond[0].index in O2 and bond[1].element.name=='hydrogen']
OH_tmp2 = [[bond[1].index,bond[0].index] for bond in top.bonds if bond[1].index in O2 and bond[0].element.name=='hydrogen']  
O2H.extend(OH_tmp1)
O2H.extend(OH_tmp2)
for i in range(len(O2H)):
   for j in range(len(O1)):
       tmp = O2H[i].copy()
       tmp.extend([O1[j]])
       tmp = np.array(tmp,dtype=int)
       triplets.append(tmp)
       OOpairs.append([O2H[i][0],O1[j]])
print('{} total OHO triplet'.format(len(triplets)))
print('Calculate OHO angles')
OOpairs = np.array(OOpairs,dtype=int)
angles = mdtraj.compute_angles(traj,triplets) * 180./np.pi # frames x triplets
angleBool = np.zeros_like(angles, dtype=bool) # all False
i = np.where(angles >= amin)
angleBool[i] = True
angleBool = np.array(angleBool, dtype=float) # turn elements to 1. or 0.
angleBool_short = []

print('Evaluate angle criterion, OHO angle must be at least {} deg'.format(amin))
# add angleBool of the angles involving same pair of oxygens
for idx,[i1,i2] in enumerate(pairs):
    Bool = np.zeros(angleBool.shape[0])
    i = np.where(OOpairs==i1)[0]
    j = np.where(OOpairs==i2)[0]
    k = np.intersect1d(i,j) # get row of triplets that involve both O1 and O2
    Bool = np.sum(angleBool[:,k],axis=1) # size = frames
    angleBool_short.append([a>0. for a in Bool]) # if at least one of the angle involving i1 and i2 oxygens satisfies the angle criterion, count it as possible H-bond
    if idx % 5 ==0:
        print("OO pair {}/{}".format(idx+1,len(pairs)), end="\r")
angleBool_short = np.array(np.array(angleBool_short).transpose(),dtype=bool) # frames x pairs
IsHbond = angleBool_short*distanceBool # satisfy both angle AND distance criteria
nHbond = np.count_nonzero(IsHbond == True,axis=1) # number of H bonds in each frame

s = 'Number of hydrogen bonds between %s oxygens and %s oxygens %5.5f +/- %5.5f'  %(atom1,atom2,np.mean(nHbond),stats.sem(nHbond))
file = open('hbond_{}.txt'.format(ext),'w')
file.write(s)
file.close()
print(s)

name='hbond_{}'.format(ext)
data = np.vstack([range(traj.n_frames), nHbond]).T
np.savetxt(name+'.dat',data,header='frame\tnHbond')

fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
axs.plot(range(traj.n_frames), nHbond, marker=None,ls='-',lw=0.75,mfc="None",ms=2)
plt.xlabel('frame')
plt.ylabel('$n_{Hbond}$')
title ='hbond {}'.format(ext)
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
