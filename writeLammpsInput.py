"""modify lammps input generated by sim to save restart and calculate rg and rdf"""

import numpy as np,argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",required=True,help="lammps input file from simt")
parser.add_argument("-s","--step",type = int, required=True, help="number of production steps")
#parser.add_argument("-np",type = int,required=True, help="number of chains")
args = parser.parse_args()

inf=args.input
outf = 'polymer.in'
thermo = False
compute = False
with open(outf,'w') as outfile:
	f=open(inf,'r')
	s=f.readline()
	while len(s):
		if '#Thermostat & time integration' in s:
			thermo = True
		elif '#restart thermo' in s:
			thermo = False
		elif '%d %d %14.7e %14.7e %14.7e %d' in s:
			compute = True
		if thermo:
			if 'thermo_style' in s:
				line='thermo_style    custom step elapsed temp ke pe etotal press vol density\n'
			elif 'thermo_modify' in s:
				line = 'restart 1000 f1.restart f2.restart'
			else:
				line = s
		elif ('thermo_style' in s or 'thermo_modify' in s) and not thermo:
				line="\n"
		elif 'dump  ' in s:
			index = s.split().index('custom')
			list = s.split()
			list[index+1] = '5000' #set stride of dump to be 5000
			index = s.split().index('id')
			list[index-1]='traj_wrapped.lammpstrj'
			list.append('\n')
			line = ' '.join(list)
	
		elif compute and len(s.split())==0:	
				line= """dump dcdtraj all dcd 5000 traj.dcd
dump_modify     dcdtraj  sort id

#compute
compute RDF_PP all rdf 100 * *
fix getRDF_PP all ave/time 5000 {} {} c_RDF_PP[*] file rdfCG_polymer_polymer.txt mode vector

#PRODUCTION RUNS \n""".format(args.step/5000,args.step)
		else:
			line = s
		outfile.write(line)
		s=f.readline()
