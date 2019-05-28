#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 21:01:12 2019

@author: nvthaomy
"""
"""Edit Lammps trajectory to use in sim"""
import sys
traj = sys.argv[1]
outTraj = traj.split('.lammpstrj')[0] + '_fixed.lammpstrj'
with open(traj,'r') as infile:
    with open (outTraj,'w') as outfile:
        removing_extra_sections = False
        line = infile.readline()
        while len(line):              
            if line.split()[0]=="ITEM:":
                if line.split()[1] == "ATOMS":
                    removing_extra_sections = True
                    id_ind = line.split().index('id') -2
                    type_ind = line.split().index('type') -2
                    x_ind = line.split().index('x') -2
                    y_ind = line.split().index('y') -2
                    z_ind = line.split().index('z') -2
                    newLine = 'ITEM: ATOMS id type x y z\n'
                    outfile.write(newLine)
                elif "TIMESTEP" in line:
                    removing_extra_sections = False    
                    outfile.write(line)
                else:
                    outfile.write(line)                   
            else:
                if removing_extra_sections:
                    cols = [line.split()[id_ind],line.split()[type_ind],line.split()[x_ind],line.split()[y_ind],line.split()[z_ind]]
                    newLine = '  '.join(cols)
                    outfile.write(newLine+'\n')
                else:
                    outfile.write(line)                            
            line = infile.readline() 