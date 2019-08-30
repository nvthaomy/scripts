#python 3
import numpy as np
import os,argparse
import MDAnalysis as mda
from MDAnalysis.topology.LAMMPSParser import DATAParser

parser = argparse.ArgumentParser()
parser.add_argument("-d", required=True,help="lammps data file")
parser.add_argument("-t",required=True, help=".dcd trajectory")
parser.add_argument("-n",required=True,type=int, help = "number of chains")
args = parser.parse_args()

LammpsData  = args.d
LammpsTrj   = args.t
NMol = args.n

u = mda.Universe(LammpsData,LammpsTrj)

#u = mda.Universe.empty()
#u.load_new()

#d = DATAParser(LammpsData)
#d.writePDB("lammps.pdb")

print (u.atoms)
list(u.atoms)

Atoms_Molecule = []
for i in range(NMol):
    molID = i+1
    Atoms_Molecule.append(u.select_atoms("resid {}".format(molID)))
#print (Atoms_Molecule)
Rg_List = []
Rg_all = [] #radius of gyration of all chain at every frame
for ts in u.trajectory:
	Rg_temp = 0	
	Rg_perframe = []
	for molecule in Atoms_Molecule:
		Rg_temp += molecule.radius_of_gyration()
		Rg_perframe.append(molecule.radius_of_gyration())
	Rg_all.append(Rg_perframe)
	Rg_List.append((u.trajectory.time, (Rg_temp/NMol)))
Rg_List = np.array(Rg_List)    
Rg_all = np.array(Rg_all)
Rg_avg = np.mean(Rg_all,axis = 0)#average Rg of each chain over all frame
print ('average Rg of each chain over all frame:')
print(Rg_avg)
import matplotlib.pyplot as plt
ax = plt.subplot(111)
ax.plot(Rg_List[:,0], Rg_List[:,1], 'r--', lw=2, label=r"$R_G$")
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)")
ax.figure.savefig("Rgyr.pdf")
plt.draw()

print ('Average Rg:')
print('average Rg over all chains, then average over frames {} +/- {}'.format(np.mean(Rg_List[:,1]),np.std(Rg_List[:,1]))) #average Rg over all chains, then average over frames
print('average Rg simultaneously over all chains and frames {} +/- {}'.format(np.mean(Rg_all),np.std(Rg_all))) #average Rg simultaneously over all chains and frames
if NMol >1:
	print('average Rg of each chain over all frames, then average over all chains {} +/- {}'.format(np.mean(Rg_avg),np.std(Rg_avg))) #average Rg of each chain over all frames, then average over all chains

with open('rg_avg.txt','w') as avg:
	avg.write('average Rg over all chains, then average over frames {} +/- {}'.format(np.mean(Rg_List[:,1]),np.std(Rg_List[:,1])))
	avg.write('\naverage Rg simultaneously over all chains and frames {} +/- {}'.format(np.mean(Rg_all),np.std(Rg_all)))
	if NMol>1:
		avg.write('\naverage Rg of each chain over all frames, then average over all chains {} +/- {}'.format(np.mean(Rg_avg),np.std(Rg_avg)))
	avg.write('\naverage Rg of each chain over all frames:')
	for i in range(NMol):
		avg.write('\n {} +/- {}'.format(Rg_avg[i],np.std(Rg_all,axis = 0)[i])) 
with open('rg.txt','w') as f:
	f.write('# frame ')
	for i in range(NMol):
		f.write('rg{} '.format(i+1))
	f.write('\n')
	for i in range(len(Rg_all)):
		rg_str = [str(x) for x in Rg_all[i]]
		f.write('{}  {}\n'.format(i+1,'  '.join(rg_str)))
