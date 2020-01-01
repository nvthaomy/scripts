import MDAnalysis as mda
from MDAnalysis.analysis import polymer
"""get persistence length of coarse-grained polymer chains assuming polymer chains have indices from 0 to np-1"""
top = 'xp0.09.data'
traj = 'xp0.09_trajwrapped.dcd'
np = 15
filter = 'type 1 2'
figname='xp0.09_persistenceL'
#-----------
u = mda.Universe(top,traj)
molecules  = u.atoms.fragments
backbones = [mol.select_atoms(filter) for mol in molecules][0:np]
lp = polymer.PersistenceLength(backbones)
lp.run()
lp.perform_fit()
LP = lp.lp
LB = lp.lb

print('persistence length: %5.5f' %LP)
print('bond length: %5.5f' %LB)
print('monomers per persistence length: %5.2f' %(LP/LB))
plot = lp.plot(figname=figname)
