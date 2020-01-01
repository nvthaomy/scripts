"""make new trajectory with head and tail monomers having a different atom type"""

import numpy as np,argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs='+',required=True,help="lammpstrj file")
parser.add_argument("-N",type = int,required=True,help="DOP")
args = parser.parse_args()
N=args.N
counter = 1
print "DOP is {}".format(N)
for txt in args.input:
    getAtomtypeIndex = True
    traj_out = txt.split('.lammpstrj')[0] + '_addedHeadTail.lammpstrj'
    with open(traj_out,'w') as outfile:
		f=open(txt,'r')
		s=f.readline()
		while len(s):
			if "ITEM: ATOMS" in s:
				readingTraj = True
                                if getAtomtypeIndex:
                                    atIndex = s.split().index('type') -2
                                    getAtomtypeIndex = False
                                    print atIndex
			elif "ITEM: TIMESTEP" in s:
				readingTraj = False
			if len(s.split()) > 2  and readingTraj and not'ITEM:' in s:
				#print "counter: {}".format(counter)
				if counter == 1 or counter == N:
					cols = s.split()
					cols[atIndex] = '2'
				        cols.append("\n")	
					s = '  '.join(cols)
				if counter== N:
					counter = 0
				counter += 1
			outfile.write(s)
			s=f.readline()

