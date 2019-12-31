import stats
import numpy as np
import matplotlib, sys, os
import matplotlib.pyplot as plt
import MDAnalysis as mda
import mdtraj as md

RgDatName = 'f0.25_RgTimeSeries'
StatOutName = 'f0.25_RgStat'
autowarmup = True
#defulat warmup samples if not using autowarmup
warmup = 100

trajFiles = ['trajectory298.dcd']
tops = ['AA12_f0.25_opc_gaff2_w0.13.parm7']
stride = 5

#names of residue in polymer chain
#resInChain = ['AHP','AP', 'ATP']
DOPs = [12]
NPs = [15]
#index of first polymer residue
res0Id = 0

ElementDictionary ={
                    "carbon": 12.01,
                    "hydrogen": 1.008,
                    "oxygen": 16.00,
                    "nitrogen": 14.001,
                    "virtual site": 1.0,
                    "virtual_site": 1.0,
                    "sodium": "na+",
                    }
#########################End of input######################
for i, trajFile in enumerate(trajFiles):
    top = tops[i]
    DOP = DOPs[i]
    NP = NPs[i]

    traj = md.load(trajFile, top=top, stride = stride)

    RgStats = []
    ReeTimeseries = [range(traj.n_frames)]
    RgTimeseries = [range(traj.n_frames)]
    header = "Frame   "
    
    StatOutExt = '_'+'.'.join(trajFile.split('.')[:-1]) + '.dat'
    txt = ""

    #get indices of residues in all chains    
    MoleculeResidueList = []
    for j in range(NP):
        resId = range(res0Id + j*DOP, res0Id + (j+1)*DOP)
        MoleculeResidueList.append(resId)

    for j,resId in enumerate(MoleculeResidueList):
        resIdLow = np.min(resId)
        resIdUp = np.max(resId)
        atom_indices = traj.topology.select('resid {} to {}'.format(resIdLow,resIdUp)) 
        mass_list = []
        for index in atom_indices:
            temp = ElementDictionary[str(traj.topology.atom(index).element)]
            mass_list.append(temp)
        mass_list = np.array(mass_list)
        Rg = md.compute_rg(traj.atom_slice(atom_indices),masses=mass_list)
 
        RgTimeseries.append(Rg.tolist())
        header += 'Rg{}   '.format(j+1)
        np.savetxt(RgDatName+StatOutExt, np.transpose(RgTimeseries), fmt = '%5.5f', header=header ) 
        
        #do stats
        file = open(RgDatName+StatOutExt,'r')
        if autowarmup:
            warmup,Data,nwarmup = stats.autoWarmupMSER(file, j+1)
            print ("Auto warmup detection with MSER-5 => ",nwarmup)
        else:
            warmup,Data = stats.extractData(file, j+1, warmup)
        (nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data, False ,False,'_{0}_mol{1}'.format(file.name,j+1))

        lines = "" 
        lines += '\n==== Rg for molecule {} ===='.format(j+1)
        lines += "\n  - Mean                    = {} +/- {}".format(mean,semcc)
        lines += "\n  - Equilibrated samples    = {}".format(nsamples)
        lines += "\n  - Correlation time        = {}".format(kappa)
        lines += "\n  - Effective # samples     = {}".format(nsamples/kappa)
        lines += "\n  - Reduced-bias variance   = {}".format(unbiasedvar)
        # note that there is no unbiased estimator for the population standard deviation. We can use sqrt(var) as a indicative estimator.
        lines += "\n  - S.D. (unbiased, biased) = {} {}".format(np.sqrt(unbiasedvar),np.std(Data,ddof=0)) # ddof is correction to 1/N...using ddof=1 returns regular reduced-bias estimator
        lines += "\n  - Min, Max                = {} {}\n".format(min,max)
        print(lines)
        txt += lines

        RgAvg = mean
        RgStd = np.sqrt(unbiasedvar)
        RgErr = semcc
        CorrTime = kappa 
        RgStats.append([RgAvg,RgStd,CorrTime,RgErr])

#        print ('The Rg for molecule {} (mean, error, std)'.format(j))
#        print ('\t{0:2.4f}\t{1:2.5f}\t{1:2.5f}'.format(RgAvg, RgErr, RgStd))

        ''' Plot the radius of gyration '''
        plt.plot(Rg, "k-")
        plt.xlabel('timestep')
        plt.ylabel('Radius-of-gryation')
        plt.savefig("Rg{}.pdf".format(j+1),bbox_inches='tight')
        plt.close()
    
    #get averages of stats
    RgStats = np.array(RgStats)
    RgAvg = np.mean(RgStats[:,0])
    RgStd = np.mean(RgStats[:,1])
    CorrTime = np.mean(RgStats[:,2])
    RgErr = np.mean(RgStats[:,3])
    RgErr_Prop = np.sqrt(np.sum(RgStats[:,3]**2))/NP
    CorrTimeErr = np.sqrt(np.var(RgStats[:,2])/len(RgStats[:,2]))
    lines = ""
    lines += '\n\n=====================\nTotal Rg average is: {0:2.3f} +/- {1:2.5f}'.format(RgAvg,RgErr)
    lines += '\nTotal Rg avg. correlation time: {0:5.4f} +/- {1:5.6f}'.format(CorrTime, CorrTimeErr)
    print(lines)
    txt += lines
    f = open(StatOutName+StatOutExt,'w')
    f.write(txt)
