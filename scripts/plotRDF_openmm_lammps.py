import matplotlib ,os
import matplotlib.pyplot as plt  
import argparse
showPlots = True
#try:
#    os.environ["DISPLAY"] #Detects if display is available
#except KeyError:
#    showPlots = False
#    matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window


"""comparing rdf from openmm and lammps, currently for simulations with lj unit
    inputs: rdf text files direcltly calculated in lammps or from mdtraj_rdf.py"""

parser = argparse.ArgumentParser() 
parser.add_argument("-i",required=True,nargs="+", help="rdf.txt files")
parser.add_argument("-p", required = True, help = "pair name") 
args = parser.parse_args() 
files = args.i
pair = args.p

plt.figure()
for i in files:
	name=i.split('.txt')
	IsLammps = False
	r=[]
	g_r=[]
	f=open(i,'r')
	line=f.readline()
	while len(line):
		if "Time-averaged" in line:
			IsLammps = True
		if not line.startswith('#'):
			if IsLammps:
				if len(line.split())>2:
					r.append(float(line.split()[1]))
					g_r.append(float(line.split()[2]))
			else:
				r.append(10*float(line.split()[0])) #"converting nm to Angstrom", to be consistent with lammps unit
				g_r.append(float(line.split()[1]))
		line=f.readline()
	plt.plot(r,g_r,label=name)

plt.xlabel('r') 
plt.ylabel('$g_r$(r)')
plt.legend(loc='best')
plt.xlim(0,10) 
plt.ylim(0)
plt.savefig('rdf-{}.png'.format(pair))
#plt.show() 

