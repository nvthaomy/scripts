import mdtraj as md
import os, sys
import argparse

import argparse as ap
parser = ap.ArgumentParser(description="mapping AA trajectory")
parser.add_argument('traj', type=str, help = "AA trajectory")
parser.add_argument('top', type=str, help = "AA topology")
parser.add_argument('ext', type=str, choices = ['lammpstrj','xyz','pdb','dcd'])
parser.add_argument('-stride', type=int, help = "stide", default = 1)
args = parser.parse_args()

traj = sys.argv[1]
top = sys.argv[2]
ext = sys.argv[3]
stride = args.stride

#outTraj = '.'.join(traj.split('.')[:-1]) +'.'+ext 
outTraj = '.'.join(traj.split('.')[:-1]) 
i=0
if os.path.exists(outTraj+'.'+ext):
    outTraj = outTraj + '_{}'.format(i)
while os.path.exists(outTraj+'.'+ext):
    i += 1
    outTraj = '_'.join(outTraj.split('_')[:-1]) + '_{}'.format(i)
outTraj = outTraj +'.'+ext

print('Loading trajectory')
traj = md.load(traj, top = top, stride = stride)
print('\nDone loading, {} frames'.format(traj.n_frames))
traj.save(outTraj)
