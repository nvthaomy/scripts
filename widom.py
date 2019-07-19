#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 14:58:16 2019

@author: nvthaomy
"""
"""Performing widom insertion for NVT"""
from itertools import islice
import numpy as np
import argparse, os
import sys
parser = argparse.ArgumentParser()
parser.add_argument("-i",required=True,help="traj file")
parser.add_argument("-at", required=True, help="atom type of inserting particle in the traj file")
parser.add_argument("-n", default = "100", help = "Number of insertions per frame")
parser.add_argument("-top", help = "topology file, e.g. pdb")
parser.add_argument("-o",default = "traj_widom", help="name of output file")
args = parser.parse_args()

def getInfo(traj):
    """get indices of atom type and x y z values and box dimensions"""
    trjFile = open(traj,'r')
    line =trjFile.readline()
    gotIndices = False
    readBox = False
    while len(line):
        if "ITEM: BOX" in line:
            readBox = True
            counter = 0
        if readBox and not "ITEM: BOX" in line: 
            if counter == 0:
                x = float(line.split()[1])-float(line.split()[0])
            elif counter == 1:
                y = float(line.split()[1])-float(line.split()[0])
            elif counter == 2:
                z = float(line.split()[1])-float(line.split()[0])
            counter += 1
            if counter > 2: 
                readBox = False        
        if "ITEM: ATOMS" in line:
            atomID_ind = (line.split()).index('id') -2
            atomtype_ind = (line.split()).index('type') -2
            x_ind = (line.split()).index('x') -2
            y_ind = (line.split()).index('y') -2
            z_ind = (line.split()).index('z') -2
            nCols = len(line.split()) -2 #number of columns
            gotIndices = True
        if gotIndices:
            line = ""
        line =trjFile.readline()
    boxL = np.array([x,y,z])
    return nCols,atomID_ind,atomtype_ind,x_ind,y_ind,z_ind,boxL

def getFrames(traj,nFrames):
    """grouping lines in traj file into a tuple"""
    trjFile = open(traj,'r')
    line =trjFile.readline()
    linesperframe = 0
    while len(line):
        if "ITEM: TIMESTEP" in line and linesperframe > 0:
            line = ""
        else:
            linesperframe += 1
            line=trjFile.readline()
    frames = [] 
    print ("Lines per frame: {}".format(linesperframe))
    with open(traj, 'r') as infile:
        lines=infile.readlines()
        for iframe in range(nFrames):
            
            frames.append(lines[iframe*linesperframe:(1+iframe)*linesperframe])
    return frames           
import mdtraj
trj = mdtraj.load(args.i, top=args.top)
nFrames = trj.n_frames
nInsertions = int(args.n)
nAtoms = trj.n_atoms
print ("Number of insertions per frame: {}".format(nInsertions))
print ("Total number of frames: {}".format(nFrames))

nCols,atomID_ind,atomtype_ind,x_ind,y_ind,z_ind,boxL = getInfo(args.i)
frames = getFrames(args.i,nFrames)
outTraj  = "{}.lammpstrj".format(args.o)
with open(outTraj,'w') as outfile:
    j = 0    
    for iframe in range(nFrames):
        linesInFrame = frames[iframe]
        linesInFrame[3] = str(nAtoms + 1)+"\n" #update number of atoms
        if np.mod(iframe,50) == 0:
                sys.stdout.write('\r')
                sys.stdout.write("Inserting test particles into frame {}".format(iframe+1))
                sys.stdout.flush()
        for i in range(nInsertions):
            linesInFrame[1] = str(j)+"\n"  #change TIMESTEP into sequential list of integer with step size =1
            randxyz = np.random.random([1,3]) * boxL
            randxyz -= boxL/2
            newline = ['1']*nCols
            newline[atomID_ind] = str(nAtoms+1)
            newline[atomtype_ind] = str(args.at)
            newline[x_ind] = str(round(randxyz[0,0],5))
            newline[y_ind] = str(round(randxyz[0,1],5))
            newline[z_ind] = str(round(randxyz[0,2],5))
            newline = " ".join(newline)+"\n"
            outfile.write("".join(linesInFrame))
            outfile.write(newline)
            j += 1
