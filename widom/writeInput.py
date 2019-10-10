import numpy as np
import argparse,os
from shutil import copyfile
"""Rewrite lammps input files from original file
 to get potential energy from a trajectory"""

parser = argparse.ArgumentParser()
parser.add_argument("-i",required=True,help="lammps input file")
parser.add_argument("-dat",required=True,help="data file")
parser.add_argument("-trj",required=True,help=".lammpstrj trajectory file")
parser.add_argument("-o",required=True,help="new input file")
parser.add_argument("-txt",required=True,help="name of output potential E file")
args = parser.parse_args()

out=args.o+".in"
cwd = os.getcwd()
input=args.i
#copyfile(os.path.join(cwd,input),os.path.join(cwd,out))

#copyfile(args.in,out)
input = open(args.i,'r')
input_str = input.read()
input_str = input_str.replace('read_data', 'read_data {} \n#read_data'.format(args.dat))
input_str = input_str.replace('\nrestart','\n#restart')
input_str = input_str.replace('\nminimize','\n#minimize')
input_str = input_str.replace('\nrun','\n#run')
input_str = input_str.replace('\ndump','\n#dump')
input_str = input_str.replace('\nwrite_data','\n#write_data')
input_str = input_str.replace('\nundump','\n#undump')
input_str = input_str.replace('\ncompute','\n#compute')
input_str = input_str.replace('\nfix getRDF','\n#fix getRDF')

#adding lines to rerun and calculate the potential energy
str = """#------Looping over trajectory with test particle and calculate potential energy-------
thermo 1
compute PE1 all pe
thermo_style custom step c_PE1
fix getPE1 all ave/time 1 1 1 c_PE1 file {}.txt format %10.10f
rerun {} every 0 dump x y z box yes replace yes""".format(args.txt,args.trj)
input_str = input_str.replace('write_data warmed_up2.data','\n{}'.format(str))
with open(out,'w') as file:
    file.write(input_str)
file.close()
