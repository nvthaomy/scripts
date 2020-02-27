import FEP_Module_short as FEP
import numpy as np
from scipy.optimize import least_squares

def GetChemPot(FEPMolNames, unit, temp, log0to1, log1to0, Plot=False):
    """get chemical potential if already have changes in energy from insertion and deletion
       and return weights"""
    FEP_Object = FEP.FEP('OB','particle')
    FEP_Object.SetTemperature(temp)
#    FEP_Object.SetNumberInsertionsPerFrame(nInsert)
#    FEP_Object.SetNumberDeletionsPerFrame(nDelete)

    dU_0to1 = np.loadtxt(log0to1, delimiter=',')[:,4] 
    dU_1to0 = np.loadtxt(log1to0, delimiter=',')[:,4]


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
    s += '\nNumber of insertion data: {}\nNumber of deletion data: {}'.format(len(dU_0to1),len(dU_1to0))
    s += "\ndF: %4.4f +/- %2.5f %s"%(dF, dF_stderr,unit)
    s += "\nBennett's Constant: %4.4f %s"%(float(dF),unit)
    print(s)
    f = open('chemicalPot.dat','w')
    f.write(s)
    f.close()

    if Plot:
        import matplotlib.pyplot as plt
        plt.bar(dU_0to1,weights_0to1,color = 'r', label='Insertion')
        plt.bar(-dU_1to0,weights_1to0,color = 'k', label='Deletion')
        plt.xlim(0.25*dF, 1.75*dF)
        plt.xlabel('dU')
        plt.ylabel('weight')
        plt.legend(loc='best')
        plt.savefig('weights.png',dpi=500,bbox_inches='tight')
        plt.show()
    return dF,dF_stderr,weights_0to1,weights_1to0

if __name__ ==  '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("MolNames",nargs='+',help="residue names to insert/delete")
    parser.add_argument("unit", type = str, help="unit of chemical potential")
    parser.add_argument("temp",type=float, help="temperature consistent with unit")
    parser.add_argument('-l0', type=str, default = 'ReRun_PotEne_State0to1_Insertions.log')
    parser.add_argument('-l1', type=str, default = 'ReRun_PotEne_State1to0_Deletions.log')
    parser.add_argument('-g', action='store_true')
    args = parser.parse_args()

    temp = args.temp
    FEPMolNames = args.MolNames
    unit = args.unit
    log0to1 = args.l0 
    log1to0 = args.l1
    dF,dF_stderr,weights_0to1,weights_1to0 = GetChemPot(FEPMolNames, unit, temp, log0to1, log1to0, Plot=args.g) 
