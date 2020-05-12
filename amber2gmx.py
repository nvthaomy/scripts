import sys
import parmed as pmd

name = sys.argv[1]
opt = int(sys.argv[2])

if opt == 0:
    print('convert amber to gromacs format')
    topExt0 = '.parm7'
    crdExt0 = '.crd'
    topExt1 = '.top'
    crdExt1 = '.gro'
    p = pmd.load_file(name+topExt0, name+crdExt0)
    p.save(name + crdExt1)
elif opt == 1:
    print('convert gromacs to amber format')
    topExt1 = '.parm7'
    crdExt1 = '.crd'
    topExt0 = '.top'
    crdExt0 = '.gro'
    p = pmd.load_file(name+topExt0)
p.save(name + topExt1)

