import mdtraj as md
import os, sys
import argparse

InTraj = os.path.abspath(sys.argv[1])
pdb = os.path.abspath(sys.argv[2])
OutTraj = sys.argv[1].split('.')[0]+'.lammpstrj'
#OutTraj2 = sys.argv[1].split('.')[0]+'.dcd'

t = md.load(InTraj,top=pdb)
t.save(OutTraj)
#t.save(OutTraj2)
