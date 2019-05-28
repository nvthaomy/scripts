#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 11:27:27 2019

@author: nvthaomy
"""
import re
import argparse
"""Convert OpenMM log.txt to a text file that can be read by stats.py"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs='+',required=True,help="log.txt file")                 
parser.add_argument("-fi",default = 'openmm', help="lammps or openmm")
parser.add_argument("-s",default = 'production',
                    help="only if input file is lammps log file: warm-up or production")
args = parser.parse_args()
def log2txt_lammps(logfiles,section,section_name): 
    start_of_section = False
    for log in logfiles:
        logOut ='log_'+ section_name
        with open(log,'r') as infile:
            with open(logOut,'w') as outfile:
                s = infile.readline()
                while s:
                    if s[0].startswith('#') and section in s:
                        start_of_section = True
                        s=infile.readline()
                    if start_of_section:
                        if s.startswith('Step'):
                            outfile.write('# '+s)
                        elif not s.startswith('Step') and not s.startswith('#') : 
                            except_num = 0
                            for value in s.split():                           
                                try: 
                                    value = float(value)                                
                                except ValueError:
                                    except_num += 1
                                    continue #trying to convert the next value in line if this value is not a float
                            if except_num == 0: #only write out lines that only contain numbers      
                                outfile.write(s)
                    if start_of_section and s.startswith('#') and ('PRODUCTION RUNS' in s or 'WARM-UP' in s):                  
                                    s = False
                    else:
                                    s=infile.readline()
                    
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
        section = 'PRODUCTION RUNS'
    elif args.s == 'warm-up':
        section = 'WARM-UP'
    log2txt_lammps(args.input,section,args.s)