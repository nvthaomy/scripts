import numpy as np
import os
from glob import glob

"""Updating residue ID so that the first resid is 1"""

files = glob('*pdb')
des = 'test'

try:
    os.mkdir(des)
except:
    pass

for fileName in files:
    FixResid = False
    SkipToNextFile = False
    file = open(fileName,'r')
    lines = file.readlines()
    linenum = 0
    newlines = []
    for line in lines:
        if 'ATOM' in line or 'HETATM' in line or 'TER' in line:
            resid = int(line[22:26])
            if linenum == 0:
                if resid == 0:
                    FixResid = True
                    newid = resid + 1
                    line = line[:22] + ' '*(4-len(str(newid))) + str(newid) + line[26:]
                    print('-- Updating residue id for {} --'.format(fileName))
                elif resid!= 0:
                    SkipToNextFile = True
            elif linenum !=0 and FixResid:
                newid = resid + 1
                line = line[:22] + ' '*(4-len(str(newid))) + str(newid) + line[26:]
            linenum += 1
        newlines.append(line)

        if SkipToNextFile:
            break

    if not SkipToNextFile:
        newfile = open(os.path.join(des,fileName.split('/')[-1]),'w')
        newfile.write(''.join(newlines))           
        newfile.close() 
