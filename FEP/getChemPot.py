import FEP_Module_short as FEP
import numpy as np
from scipy.optimize import least_squares

"""get chemical potential if already have changes in energy from insertion and deletion
   and return weights"""

FEPMolNames = ['Na+', 'Cl-']
nInsert = 20
nDelete = 5
nFrames = 950
unit = 'kJ/mol'

FEP_Object = FEP.FEP('OB','particle')
FEP_Object.SetTemperature(298.)
FEP_Object.SetNumberInsertionsPerFrame(nInsert)
FEP_Object.SetNumberDeletionsPerFrame(nDelete)

dU_0to1 = np.loadtxt('ReRun_PotEne_State0to1_Insertions.log',delimiter=',')[:,4] 
dU_1to0 = np.loadtxt('ReRun_PotEne_State1to0_Deletions.log',delimiter=',')[:,4]


#---------#
beta = 1./FEP_Object.ThermalEnergy # inverse KJ/mole
print('beta {}'.format(beta))

ExpdU_0to1 = np.exp(-1*beta*dU_0to1/2.)
ExpdU_1to0 = np.exp(-1*beta*dU_1to0/2.)

opt = least_squares(FEP_Object.BennettsCostFunction,FEP_Object.BennettsConstant,args=(dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta))

OptimalBennettsConstant = opt.x
dF,dF_stderr,weights_0to1,weights_1to0 = FEP_Object.BennettsDeltaFreeEnergy(OptimalBennettsConstant,dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta)

FEP_Object.BennettsConstant = float(OptimalBennettsConstant)


A1 = np.column_stack((dU_0to1,weights_0to1))
A2 = np.column_stack((dU_1to0,weights_1to0))
np.savetxt('dU0to1.dat',A1, header='dU weight')
np.savetxt('dU1to0.dat',A2, header='dU weight')

s = 'Excess chemical potential calculation for {}'.format(FEPMolNames)
s += '\nNumber of insertions: {}\nNumber of deletions: {}\nNumber of traj frames: {}'.format(nInsert,nDelete,nFrames)
s += "\ndF: %4.4f +/- %2.5f %s"%(dF, dF_stderr,unit)
s += "\nBennett's Constant: %4.4f %s"%(float(dF),unit)
print(s)
f = open('chemicalPot.dat','w')
f.write(s)

