#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 11:27:27 2019

@author: nvthaomy
"""
import re
import argparse

"""Convert OpenMM log.txt and lammps trj (exported by sim) to a text file that can be read by stats.py"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs='+',required=True,help="log.txt file")                 
parser.add_argument("-fi",default = 'openmm', help="lammps or openmm")
parser.add_argument("-s",default = 'production',
                    help="only if input file is lammps log file: warm-up or production")
args = parser.parse_args()
def log2txt_lammps(logfiles,section,section_name): 
    start_of_section = False
    for log in logfiles:
        step =[]
        TotEng=[]
        KinEng=[]
        PotEng=[]
        Temp=[]
        E_bond=[]
        E_angle=[]
        E_dihed=[]
        E_vdwl=[]
        logOut ='log_'+ section_name +'.txt'
        with open(log,'r') as infile:
            s = infile.readline()
            while s:
		#if s[0].split().startswith('#') and section in s:
                if '#' in s  and section in s:
                    start_of_section = True
                    s=infile.readline()
                if start_of_section:
                    if '---------------- Step' in s:
                        step.append(s.split()[2])
                    elif 'TotEng' in s:
                        TotEng.append(s.split()[2])
                        KinEng.append(s.split()[5])
                        PotEng.append(s.split()[8])
                    elif 'Temp' in s:
                        Temp.append(s.split()[2])
                        E_bond.append(s.split()[5])
                        E_angle.append(s.split()[8])
                    elif 'E_dihed' in s:
                        E_dihed.append(s.split()[2])
                        E_vdwl.append(s.split()[5])
                if start_of_section and s.startswith('#') and ('run production' in s or 'run equilibration' in s):                  
                                s = False
                else:
                                s=infile.readline()                 
        with open(logOut,'w') as outfile:
            outfile.write('# Step TotEng KinEng PotEng Temp E_bond E_angle E_dihed E_vdwl\n')
            for ind,val in enumerate(step):
                outfile.write('{} {} {} {} {} {} {} {} {}\n'.format(step[ind],
                              TotEng[ind],KinEng[ind], PotEng[ind], Temp[ind], 
                              E_bond[ind], E_angle[ind], E_dihed[ind], E_vdwl[ind]))
def log2txt_openmm(logfiles):
    """reading multiple log files and convert to readable text files"""
    for log in logfiles:
        logOut = 'data'+log[len('log'):]
        with open(log,'r') as infile:
            with open(logOut,'w') as outfile:
                s = infile.readline()
                while s:
                    if s[0].startswith('#'):
                        cols = [x for x in re.split('"|\t|\n',s) if x != ''] #removing " and white space between " "
                        new_cols = []
                        for i in cols:
                            i = i.replace(' ','_')
                            new_cols.append(i)
                        if 'Time_Remaining' in new_cols:
                            time_on = "True"
                            index_time = new_cols.index('Time_Remaining')
                            del new_cols[index_time]
                        s = '   '.join(new_cols)
                        outfile.write(s)
                        outfile.write('\n')
                    else:
                        s = s.replace('%','')
                        if time_on:
                            cols = s.split()
                            del cols[index_time-1]
                            s = '  '.join(cols)
                        outfile.write(s)
                        outfile.write('\n')
                    s = infile.readline()
if args.fi == 'openmm':
    log2txt_openmm(args.input)
elif args.fi == 'lammps':
    if args.s == 'production':
        section = 'run production'
    elif args.s == 'warm-up':
        section = 'run equilibration'
    log2txt_lammps(args.input,section,args.s)
