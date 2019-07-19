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
parser.add_argument("-nd", help = "Turn on nondimensionalization", action="store_true")
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
        time_on = False
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
def log2txt_openmm_nondim(logfiles,epsilon,sigma,mass):
    """read multiple log files, nondimensionalize and convert to readable text files"""
    for log in logfiles:
        N_av = 6.022140857*10**23 #/mol
        kb = 1.380649*10**(-23)* N_av * 10**-3 #KJ/kelvin/mol
        logOut = 'data'+log[len('log'):]
        time_on = False
        observables = {'PE':'Potential_Energy','TotE':'Total_Energy',
                       'T':'Temperature','V':'Box_Volume','rho':'Density',
                       'KE':'Kinetic_Energy'}
        observablesInd = {'PE':[],'TotE':[],'KE':[],'T':[],'V':[], 'rho':[]}
        with open(log,'r') as infile:
            with open(logOut,'w') as outfile:
                s = infile.readline()
                while s:
                    if s[0].startswith('#'):
                        cols = [x for x in re.split('"|\t|\n',s) if x != ''] #removing " and white space between " "
                        new_cols = []
                        for i in cols:
                            if len(i.split())>1:
                                i = i.split()[:-1] #obmitting the unit
                                i = '_'.join(i)
                            new_cols.append(i)
                        if 'Time_Remaining' in new_cols:
                            time_on = "True"
                            index_time = new_cols.index('Time_Remaining')
                            del new_cols[index_time]
                        for val in observables:
                            if observables[val] in new_cols:
                                 observablesInd[val] = new_cols.index(observables[val]) -1 # get index of observable, -1 to account for  the # symbol
                        s = '   '.join(new_cols)
                        outfile.write(s)
                        outfile.write('\n')
                    else:
                        s = s.replace('%','')
                        cols = s.split()
                        if time_on:                 
                            del cols[index_time-1]
                            s = '  '.join(cols)
                        #convert to dimensionless unit:
                        for val in observables:
                            if observablesInd[val]:
                                ind = observablesInd[val]
                                dim_val = float(cols[ind])
                                if val in ['PE','TotE','KE']:
                                    nondim_val = dim_val/epsilon #kj/mol / kj/mol
                                elif val == 'T':
                                    nondim_val = dim_val * kb / epsilon
                                elif val == 'V':
                                    nondim_val = dim_val/sigma #nm/nm
                                elif val == 'rho':
                                    nondim_val = dim_val * sigma**3 *(10**-7)**3 /mass * N_av #particles
                                cols[ind] = str(round(nondim_val,4))
                        s = '  '.join(cols)
                        outfile.write(s)
                        outfile.write('\n')
                    s = infile.readline()    
if args.fi == 'openmm':   
    if args.nd:
        epsilon = float(raw_input('\nEnter energy scale (kJ/mol) :'))
        sigma = float(raw_input('\nEnter length scale (nm) :'))
        mass = float(raw_input('\nEnter mass scale (g/mol) :'))
        log2txt_openmm_nondim(args.input,epsilon,sigma,mass)
    else:
        log2txt_openmm(args.input)
elif args.fi == 'lammps':
    if args.s == 'production':
        section = 'PRODUCTION RUNS'
    elif args.s == 'warm-up':
        section = 'WARM-UP'
    log2txt_lammps(args.input,section,args.s)
