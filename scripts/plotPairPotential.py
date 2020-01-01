"""Plot pair potential from lammps pair file
and calculate the excluded volume and the energy integral"""
import numpy as np
import argparse, os, sys
import matplotlib.pyplot as plt
import matplotlib
from scipy import integrate
showPlots = True
try:
    os.environ["DISPLAY"] #Detects if display is available
except KeyError:
    showPlots = False
    matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window

parser = argparse.ArgumentParser()
parser.add_argument("-i", nargs = "+", required = True, help = "lammps pair files")
parser.add_argument("-l", nargs = "+", required = True, help = "legends for plotting")
parser.add_argument("-n", required = True, help = "name of data set, use in plot name")
parser.add_argument("-x", nargs = "+", help = "x values for excluded volume and energy integral plots, same length as number of input files")
parser.add_argument("-xl", help = "x label for excluded volume and energy integral plots, same length as number of input files")
args = parser.parse_args()
#  python plotPairPotential.py -i xp0.024_wrapped_CGMap_3_Spline_NoP.pair xp0.024_wrapped_CGMap_4_Spline_NoP.pair xp0.024_wrapped_CGMap_6_Spline_NoP.pair xp0.024_wrapped_CGMap_9_Spline_30knots_cut15_NoP.pair xp0.024_wrapped_CGMap_12_Spline_NoP.pair  -l 3:1 4:1 6:1 9:1 12:1 -n xp0.024 -x 3 4 6 9 12 -xl "Mapping Ratio"

command_args = "-i "
for file in args.i:
    command_args += "{} ".format(file)
command_args += "-l "
for l in args.l:
    command_args += "{} ".format(l)
command_args += "-x "
for x in args.x:
    command_args += "{} ".format(x)
command_args += "-xl \"{}\" ".format(args.xl)    
command_args += "-n \"{}\" ".format(args.n)

files = args.i
if not args.x:
    xs = range(len(files))
    xlabel = 'File index'   
else: 
    xs = args.x
    xlabel = args.xl 
#make empty arrays (will have dimension of #files by Ntable) of distance, potential, force
rs = []
us = []
fs = []
Ntable = [] #number of data points in each file

#make empty arrays of excluded volume and energy integral
vs = []
es = []


rmax = 0
umax = 0
umin = 0

for file in files:
    print (file)
    infile = open(file,'r')
    s = infile.readline()
    r = []
    u = []
    f = []
    while len(s):
        #terminate when reading the NonBondNull section
        if "NonBondNull" in s:
            break
        #start to import when the first character is a float
        try:
            float(s.split()[0])
            vals = s.split()
            r.append(float(vals[1]))
            u.append(float(vals[2]))
            f.append(float(vals[3]))
            s = infile.readline()
        except:
            s = infile.readline()
    rs.append(r)
    us.append(u)
    fs.append(f)
    Ntable.append(len(r))
    #determine ranges for plotting
    if max(r) > rmax:
        rmax = max(r)
    if max(u) > umax:
        umax = max(u)
    if min(u) < umin:
        umin  = min(u)

#calculating the excluded volume
for i,u in enumerate(us):
    u = np.array(u,float)
    r = np.array(rs[i])
    mayerf =  np.exp(-u) - 1
    r2 = r**2
    y =  4 * np.pi * r2 * (-mayerf)
    v = integrate.simps(y, r)
    y = 4 * np.pi * r2 * u
    e = integrate.simps(y, r)
    vs.append(v)
    es.append(e)

#write to file
outfile = open('{}_plotPairPotential.txt'.format(args.n),'w')
s = "#{} excludedV Eintegral\n".format(xlabel)
for i, x in enumerate(xs):
    s += '{x} {v} {e}\n'.format(x=x,v=vs[i],e=es[i])
s += "#ran with arguments:\n#{}".format(command_args)
outfile.write(s)
outfile.close()
#plot pair potential
plt.figure(figsize=[3,2])
for i,r in enumerate(rs): 
    label = args.l[i]
    plt.plot(r,us[i], lw = 1., label = label)
plt.ylabel('Pair potential')
plt.xlabel('r')
plt.xlim(0,rmax)
plt.ylim(umin * 1.1, umax *1.1)
plt.legend(loc='best')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
plt.rc('legend', fontsize=7)    
plt.savefig('{}_PairPotential1.png'.format(args.n),dpi=500,bbox_inches="tight")
#plt.show()

plt.ylim(umin * 1.1, 20)
plt.savefig('{}_PairPotential2.png'.format(args.n),dpi=500,bbox_inches="tight")
plt.show()

#plot excluded volume
plt.figure(figsize=[3,2])
try:
    xs = [float(x) for x in xs]
    plt.plot(xs,vs, 'bo:',markersize = 3, lw = .75, mfc = "None")  
    plt.plot(xs,np.zeros(len(xs)),'k-', lw = 0.5)
except:
    plt.plot(range(len(xs)),vs, 'bo:',markersize = 3, lw = .75, mfc = "None")
    plt.plot(range(len(xs)),np.zeros(len(xs)),'k-', lw = 0.5)
    plt.xticks(np.arange(len(xs)),xs)
plt.ylabel('Excluded Volume')
plt.xlabel(xlabel)
plt.ylim(min(vs)-0.2*abs(min(vs)),max((max(vs)+0.2*abs(max(vs))),20))
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
plt.rc('legend', fontsize=7)
plt.savefig('{}_ExcludedV.png'.format(args.n),dpi=500,bbox_inches="tight")
plt.show()

#plot energy integral 
plt.figure(figsize=[3,2])
try:
    xs = [float(x) for x in xs]
    plt.plot(xs,es, 'ro:',markersize = 3, lw = .75, mfc = "None")
    plt.plot(xs,np.zeros(len(xs)),'k-', lw = 0.5)
except:
    plt.plot(range(len(xs)),es, 'ro:',markersize = 3, lw = .75, mfc = "None")
    plt.plot(range(len(xs)),np.zeros(len(xs)),'k-', lw = 0.5)
    plt.xticks(np.arange(len(xs)),xs)
plt.ylabel('Energy Integral')
plt.xlabel(xlabel)
plt.ylim(min(es)-0.2*abs(min(es)),max((max(es)+0.2*abs(max(es))),20))
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
plt.rc('legend', fontsize=7)
plt.savefig('{}_EnergyIntegral.png'.format(args.n),dpi=500,bbox_inches="tight")
plt.show()
