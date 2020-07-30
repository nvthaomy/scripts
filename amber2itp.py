import sys
import parmed as pmd

name = sys.argv[1]
print('convert amber to gromacs .top and .itp')
topExt0 = '.parm7'
crdExt0 = '.crd'
p = pmd.load_file(name+topExt0, name+crdExt0)
p.save(name + '.top',parameters=name+'.itp',combine=None)

