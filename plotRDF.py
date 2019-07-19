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
parser.add_argument("-i",required=True,help="rdf data file from lammps")
parser.add_argument("-p", required = True, help = "pair name")
parser.add_argument("-N", default = 200, help = "number of bins", type=int)
args = parser.parse_args()
Nbins = args.N
rdfFile = args.i
pair = args.p
f = open(rdfFile,'r')
line = f.readline()
r = np.zeros(Nbins)
g_r = np.zeros(Nbins)
coord = np.zeros(Nbins)
getInfo = False
counter = 0
while len(line):
	if '# Row' in line:
		getInfo = True
	if getInfo and not line.startswith("#") and len(line.split())>2:
		ind = int(line.split()[0])-1
		r[ind] = float(line.split()[1])
		g_r[ind] += float(line.split()[2])
		coord[ind] += float(line.split()[3])
		if ind == 0:
			counter += 1
	line = f.readline()

print "Number of frames to average: {}".format(counter)
g_r = g_r/counter
coord = coord/counter
plt.figure()
plt.plot(r,g_r)
plt.xlabel('r')
plt.ylabel('g_{}(r)'.format(pair)) 
plt.title('g_{}(r)'.format(pair),loc='center') 
plt.ylim(0)
plt.xlim(0)
plt.savefig('g_{}.png'.format(pair),dpi=500,transparent=True)
plt.show()

plt.figure()
plt.plot(r,coord)
plt.xlabel('r')
plt.ylabel('coordination_{}(r)'.format(pair)) 
plt.ylim(0)
plt.xlim(0)
plt.title('Coordination_{}(r)'.format(pair),loc='center') 
plt.savefig('coordination_{}.png'.format(pair),dpi=500,transparent=True)
plt.show()  
