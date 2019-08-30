import mdtraj as md
import numpy,argparse
"""slice trajectory by atom indices and calculate the radius of gyration for each polymer chain
   assuming the first N*np atoms are polymers and the remainders are solvent molecules"""
parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True,help="trajectory")
parser.add_argument("-N",type = int,required=True,help="DOP")
parser.add_argument("-np",type = int,required=True,help="number of chains")
parser.add_argument("-top",help="topology or pdb file")
args = parser.parse_args()

trj =args.i
pdb = args.top
N = args.N
np = args.np
if trj.split('.')[1] != 'pdb':
	t=md.load(trj,top=pdb)
else:
	t=md.load(trj)

d={} #library of trajectories with single chain
rg = {} #library of rg for every chain in every frame
rg_avg = [] #list of rg for every chain averaged over all frames
std = []

id_thusfar = 0
with open('rg_mdtraj.txt','w') as f:
        f.write('# rg std')
        for index,i in enumerate(range(np)):
                ind=range(id_thusfar,id_thusfar+N)
                key="poly{}".format(index)
                d.update({key:t.atom_slice(ind)})
                rg.update({"rg{}".format(index):10*md.compute_rg(d[key])}) #multiply by 10 when reading from lammps trajectory in lj unit
                rg_avg.append(numpy.mean(rg["rg{}".format(index)]))
                std.append(numpy.std(rg["rg{}".format(index)]))
                id_thusfar += N
                f.write('\n{}  {}'.format(rg_avg[index],std[index]))
        f.close()
print (rg)
