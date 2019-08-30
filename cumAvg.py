#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:21:22 2019

@author: nvthaomy
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import argparse

showPlots = True
try:
    os.environ["DISPLAY"] #Detects if display is available
except KeyError:
    showPlots = False
    matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window

parser = argparse.ArgumentParser()
parser.add_argument("-f",required=True,help="data file")
parser.add_argument("-o", nargs = '+', type = str, required=True,help="observables")
args = parser.parse_args()

f = open(args.f,'r')
s = f.readline()
o = []
for i in args.o:
    o.append([]) #empty list of observable values
    

while len(s):
    if s.startswith('#'):
                    obs_indices = []
                    for obs in args.o:
                        obs_indices.append(s.split().index(obs) - 1)
    else:
        cols = s.split()
        for i, obs in enumerate(args.o):
            column_num = obs_indices[i]
            o[i].append(float(cols[column_num]))
    s = f.readline()
f.close()
#getting cumulative avg:
o_avg = np.array(o)
for i,obs in enumerate(args.o):
    o_avg[i] = np.cumsum(o[i])
    n = np.array(range(1,len(o_avg[i])+1))
    o_avg[i] = o_avg[i]/n
    
#making subplots
nrows= 2
ncols = int(np.ceil(len(args.o)/nrows))
plt.figure()
for i,obs in enumerate(args.o):
    plt.subplot(nrows,ncols,i+1)
    y = o_avg[i]
    x = range(0,len(y))
    plt.plot(x,y)
    plt.title('Cumulative average of {}'.format(obs))
    plt.xlabel("Sample index")
    plt.ylim(0.99*np.min(y[100:]),1.1*np.max(y[100:]))
string = "_".join(args.o)
plt.subplots_adjust(hspace=1)
plt.savefig("CumAvg_{}.png".format(string),dpi=500)
plt.show()

                    
                        

