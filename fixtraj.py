"""make new trajectory with head and tail monomers having a different atom type"""

import numpy as np,argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs='+',required=True,help="lammpstrj file")
parser.add_argument("-N",type = int,required=True,help="DOP")
args = parser.parse_args()
N=args.N
counter = 1
print "DOP is {}".format(N)
with open('traj_new.lammpstrj','w') as outfile:
	for txt in args.input:
		f=open(txt,'r')
		s=f.readline()
		while len(s):
			if "ITEM: ATOMS" in s:
				readingTraj = True
			elif "ITEM: TIMESTEP" in s:
				readingTraj = False
			if len(s.split()) == 6 and readingTraj:
				#print "counter: {}".format(counter)
				if counter == 1 or counter == N:
					cols = s.split()
					cols[5] = '2'
				        cols.append("\n")	
					s = '  '.join(cols)
				if counter== N:
					counter = 0
				counter += 1
			outfile.write(s)
			s=f.readline()

