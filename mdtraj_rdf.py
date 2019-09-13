import mdtraj,matplotlib ,os
import matplotlib.pyplot as plt  
import os
import argparse
"""get polymer-polymer and polymer-solvent RDF using mdtraj
   currently assuming """
parser = argparse.ArgumentParser()
parser.add_argument("-i",required=True,help="trajectory in pdb format")
parser.add_argument("-N", required = True, help = "DOP")
parser.add_argument("-np", required = True, help = "number of chains")
parser.add_argument("-L",type=float,help = "box length")
args = parser.parse_args()

traj = args.i
pair1='polymer-solvent'
pair2='polymer-polymer'
L = args.L
f=open(traj,'r')
s=f.readline()
line_num = 1
Box = False 
while len(s):
	if 'CRYST1' in s:
		Box = True
		break
	if line_num > 10 and not Box:
		f.close()
		new_t = 'trajectory_box.pdb'
		new_traj = open(new_t,'w')
		f=open(traj,'r')
		line=f.readline()
		line_old = ''
		while len(line):
			if "REMARK" in line_old: 
				new_traj.write('CRYST1   {:.3f}   {:.3f}   {:.3f}  90.00  90.00  90.00 P 1           1\n'.format(L,L,L))
			new_traj.write(line)
			line_old = line
			line=f.readline()
		traj = new_t
	s=f.readline()
	line_num += 1
t=mdtraj.load(traj)
#t.save('trajectory.xyz')
top = t.topology
P=[]
S=[]
for atom_index in range(t.n_atoms):
	atom_name = str(top.atom(atom_index))
	if atom_name =="Pol1-P":
		P.append(atom_index)
	elif atom_name =="Sol1-S":
		S.append(atom_index)
for pair in [pair1,pair2]:
	if pair == 'polymer-solvent':
		pairs = top.select_pairs(selection1=P, selection2=S)
	elif pair == 'polymer-polymer':
		pairs = top.select_pairs(selection1=P, selection2=P)
	r,g_r=mdtraj.compute_rdf(t,pairs=pairs,r_range=(0,1),n_bins=200)
	name='rdf_{}.txt'.format(pair)
	f=open(name,'w')
	f.write('# r g_r_{}\n'.format(pair))
	for i in range(len(r)):
		f.write('{}    {}\n'.format(r[i],g_r[i]))
	f.close()

#plt.figure()
#plt.plot(r,g_r,label=pair1) 

#plt.plot(r,g_r,label=pair2)
#plt.xlabel('r') 
#plt.ylabel('g_{}(r)')
#plt.xlim(0) 
#plt.ylim(0)
#plt.show() 
