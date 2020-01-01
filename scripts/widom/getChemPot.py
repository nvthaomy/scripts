#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:12:12 2019

@author: nvthaomy
"""
import os
import numpy as np
import argparse
"""Calculating excess chemical potential from inserting a test particle
   Input: two data files from lammps, one with potential energy before insertion
           and one with potential energy after insertion"""
           
parser = argparse.ArgumentParser()
parser.add_argument("-u0",required=True,help="data file with u0 values")
parser.add_argument("-u1",required=True,help="data file with u1 values")
parser.add_argument("-n",required=True, type=int,help="Number of insertions per frame")
args = parser.parse_args()

nInsertions = args.n
f = open(args.u0,'r')
u0 = []
line=f.readline()
while len(line):
    if not line.startswith('#'):
                       u0.append(float(line.split()[1]))
    line = f.readline()       
nFrames = len(u0)

f=open(args.u1,'r')
u1 = []
line = f.readline()
u1_perframe = []
while len(line):       
        if not line.startswith('#'):
                               u1_perframe.append(float(line.split()[1]))
        if len(u1_perframe) == nInsertions:
            u1.append(u1_perframe)
            u1_perframe = []
        line = f.readline()
u0 = np.array(u0)
u0 = np.reshape(u0,[len(u0),1])
u1 = np.array(u1)
A = np.mean(np.exp((-u1+u0)),1) #average of exp(-deltaU*) over number of insertions
chemPot = -np.log(A)
data = np.hstack( [ chemPot[:,None], u0, u1 ] )
np.savetxt("insertion.dat", data, header="ChemPot\t\tU0\t\tU1" )
B = np.mean(A) #average of exp(-deltaU*) over traj frames

dF= -np.log(B) 

fname = 'widom_'
fId = 0
fext = 'Ninsertion_{}'.format(str(int(args.n)))
filename = fname + fext + str(fId) +'.txt'
while os.path.isfile(filename) == True:
    fId += 1
    filename = fname + fext + str(fId) +'.txt'
print ('Excess chemical potential from {} insertions: {}'.format(args.n, dF))
f = open(filename,'w')
f.write('Excess chemical potential from {} insertions: {}'.format(args.n, dF))                       
