import mdtraj as md
import os, sys
import argparse

InTraj = os.path.abspath(sys.argv[1])
OutTraj = sys.argv[1].split('.')[0]+'.xyz'

t = md.load(InTraj)
t.save(OutTraj)

