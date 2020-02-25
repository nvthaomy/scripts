
import mdtraj as md
import os, sys
import argparse

import argparse as ap
parser = ap.ArgumentParser(description="mapping AA trajectory")
parser.add_argument('top', type=str, help = "AA topology")
parser.add_argument('ext', type=str, choices = ['lammpstrj','xyz','pdb','dcd'])
parser.add_argument('-traj', type=str, nargs='+',  help = "AA trajectory")
parser.add_argument('-stride', type=int, help = "stide", default = 1)
args = parser.parse_args()

trajs = args.traj
top = sys.argv[1]
ext = sys.argv[2]
stride = args.stride
cwd = os. getcwd()

outTraj = './traj_combined.'+ext 
trajDict = {}
trajAll = 0
print('Loading trajectory')
for i,traj in enumerate(trajs):
    t = md.load(traj, top = top, stride = stride)
    trajDict.update({'t{}'.format(i):t})
    if trajAll == 0:
        trajAll = t
    else:
        trajAll += t
print('\nDone loading and combining trajectories, Total frames: {}'.format(trajAll.n_frames))
trajAll.save(outTraj)
