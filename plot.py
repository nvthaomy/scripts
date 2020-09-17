from scipy.interpolate import interp1d
import os, sys, re
import numpy as np
import mdtraj as md
import matplotlib
sys.path.append('/home/mnguyen/bin/scripts/')
import stats

showPlots = True
try:
  os.environ["DISPLAY"] #Detects if display is available
except KeyError:
  showPlots = False
  matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
import matplotlib.pyplot as plt
#plt.style.use('seaborn-dark')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
colors = ['#6495ED','r','#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#949c2d', '#a34a17','#c43b99','#949c2d','#1E90FF']

######

import argparse as ap
parser = ap.ArgumentParser()
parser.add_argument('-c', required=True, type=int, nargs = '+', help='column indices for x and y values')
parser.add_argument('-f', required=True, nargs ='+', help='data files, plot as series on the same graph if providing more than 1 data files')
parser.add_argument('-xmax', type=float, help='x max')
parser.add_argument('-xmin', type=float, help='x min')
parser.add_argument('-ymax', type=float, help='y max')
parser.add_argument('-ymin', type=float, help='y min')
parser.add_argument('-x', type=str, default = 'x', help='x label')
parser.add_argument('-y', type=str, default = 'y', help='y label')
parser.add_argument('-legend', type = str, default = None, nargs = '+', help = 'legends')
parser.add_argument('-n', type=str, default = None, help='plot name')
parser.add_argument('-m', type=str, nargs ='+', default = 'None', help='marker style')
parser.add_argument('-ls', type=str, nargs ='+', default = '-', help='line style')

args = parser.parse_args()

xlabel = args.x
ylabel = args.y
datafiles = args.f
nseries = len(args.f)
legends = args.legend

if legends == None:
    legends = ['Series {}'.format(a) for a in range(1,nseries+1)]
else:
    if not isinstance(legends,list):
        legends = [legends]
if len(args.c) != 2*nseries:
    raise Exception ('Number of columns, -c, must be 2 * number of series to plot')
if len(legends) != nseries:
    raise Exception ('Number of legends, -legend, must be the same as the number of series to plot')

plotName = args.n
marker = args.m
ls = args.ls
if not isinstance(marker,list):
    marker = [marker]
if not isinstance(ls, list):
    ls = [ls]
if len(marker) < nseries:
    marker = [marker[0]]*nseries
if len(ls) < nseries:
   ls = [ls[0]]*nseries

xs = []
ys = []
for i,datafile in enumerate(datafiles):
    if nseries == 1:
        if plotName == None:
            plotName = datafile + '_col_{}_{}'.format(args.c[0],args.c[1])
    else:
        if plotName == None:
            plotName = 'plot'
    x = np.loadtxt(datafile)[:,args.c[0]]
    y = np.loadtxt(datafile)[:,args.c[1]]
    xs.append(x)
    ys.append(y)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
for i,x in enumerate(xs):
    y = ys[i]
    if nseries == 1:
        ax.plot(x,y, marker=marker[i], ms=4,ls=ls[i],lw=1)
    elif nseries > 1:
        ax.plot(x,y, marker=marker[i], ms=4,ls=ls[i],lw=1,label = legends[i])
plt.xlabel(xlabel)
plt.ylabel(ylabel)
if nseries > 1:
    plt.legend(loc='best',prop={'size':5})
if args.xmin != None and args.xmax != None:
    plt.xlim(args.xmin,args.xmax)
else:
    if args.xmin:
        plt.xlim(xmin=args.xmin)
    elif args.xmax:
       plt.xlim(xmax=args.xmax)
if args.ymin != None and args.ymax != None:
    plt.ylim(args.ymin,args.ymax)
else:
    if args.ymin:
        plt.ylim(ymin=args.ymin)
    elif args.ymax:
       plt.ylim(ymax=args.ymax)
title = plotName
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
