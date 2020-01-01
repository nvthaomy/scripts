import glob
import sys,os
import mdtraj as md
"""1st argument: trajectory file,
   2nd argument: first frame to include in the final trajectory
   Convert .dcd trajectory to lammps traj"""
parm7 = glob.glob('*.parm7')[0]
trj = sys.argv[1]
i = int(sys.argv[2]) #first frame to include in the final trajectory
trajOut=trj.split(trj[trj.index('.'):])[0]+'.lammpstrj'
trajOutnc=trj.split(trj[trj.index('.'):])[0]+'.nc'
t=md.load(trj,top = parm7)
print ('Converting the last {} frames of the trajectory'.format(len(t)-i))
#t[i:].save(trajOut)
t[i:].save(trajOutnc)
