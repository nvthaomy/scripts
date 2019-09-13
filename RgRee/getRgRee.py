"""Calculate Rg, Ree from a trajectory"""

CalcRg          		= True
CalcRee         		= True
CnvDATA2PDB				= True	# Convert lammps.data to .pdb using MDAnalysis
Make_molecules_whole    = True  # Uses MDTraj to make molecules whole
SaveCalculationsToFile 	= True

import numpy as np
import matplotlib, sys, os
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import MDAnalysis as mda
import mdtraj as md
from pymbar import timeseries
from HistogramTools import HistogramRee
import stats_TrajParse as stats_mod
import stats

#=========================INPUT===========================
data_files = ['np01.data','np04.data','np08.data','np15.data','np40.data','np50.data','np65.data','np80.data'] #lammps data file
traj_files = ['AAtraj_np01_mapped.lammpstrj','AAtraj_np04_mapped.lammpstrj','AAtraj_np08_mapped.lammpstrj','AAtraj_np15_mapped.lammpstrj','AAtraj_np40_mapped.lammpstrj','AAtraj_np50_mapped.lammpstrj','AAtraj_np65_mapped.lammpstrj','AAtraj_np80_mapped.lammpstrj'] # lammpstrj or dcd
natoms_list = len(data_files)*[8] #number of atoms in each polymer
#----------------------------------------------------------
for index, traj_file in enumerate(traj_files):
    SysName = ''.join(traj_file.split('.')[:-1])
    sys.stdout.write('Calculating Rg and Ree for {}'.format(SysName))

    def statistics_py(DataFilename,Col):
        ''' Do data using stats.py '''
        # Calculate using stats.py
        f = open(DataFilename,"rw")
        warmupdata, proddata, idx = stats.autoWarmupMSER(f,Col)
        nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor = stats.doStats(warmupdata,proddata)

        return ([nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor])
        
    def pymbar_statistics(dataset):
        ''' Do PyMbar Analysis on the data using timeseries '''
        dataset = np.asarray(dataset).flatten()
        dataset_temp = dataset
        pymbar_timeseries = (timeseries.detectEquilibration(dataset)) # outputs [t0, g, Neff_max] 
        t0 = pymbar_timeseries[0] # the equilibrium starting indices
        Data_equilibrium = dataset#[t0:]
        g = pymbar_timeseries[1] # the statistical inefficiency, like correlation time
        indices = timeseries.subsampleCorrelatedData(Data_equilibrium, g=g) # get the indices of decorrelated data
        dataset = Data_equilibrium[indices] # Decorrelated Equilibrated data
        dataset = np.asarray(dataset)
        P_autocorrelation = timeseries.normalizedFluctuationCorrelationFunction(dataset, dataset, N_max=None, norm=True)
        # Calculate using stats.py
        warmupdata = 0
        ''' Use stats.py to calculate averages '''
        nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor = stats_mod.doStats(warmupdata,dataset_temp)

        return ([np.mean(dataset),np.var(dataset),np.sqrt(np.divide(np.var(dataset),len(dataset))),len(indices), len(Data_equilibrium), g, mean, semcc, kappa])
    
    import stats_TrajParse as stats
    # Walk through and find trajectory for specific system
    LammpsData  = data_files[index]
    LammpsTrj = traj_file
    natoms  = natoms_list[index] #number of atoms in a polymer
    DOP = natoms #DOP is natoms for coarse-grained/model system
    SaveFilename = SysName+'_RgReeData'
    if SaveCalculationsToFile: os.mkdir(SaveFilename)

    # Converts LAMMPS.data to .pdb with structure information
    u = mda.Universe(LammpsData)
    gr = u.atoms
    gr.write(SysName+'.pdb')
    top_file = SysName+'.pdb'
    
    ''' Load in trajectory file '''
    traj = md.load(LammpsTrj,top=top_file)#, atom_indices=atoms_list)
    print ("Multiply coordinates and unit cell lengths by 10 to convert to lj unit")
    traj.xyz *= 10 #multiply by 10 to convert to lj unit
    traj.unitcell_lengths *= 10 #multiply by 10 to convert to lj unit 

    resID = []
    for  res_id,residue in enumerate(traj.topology.residues):
        if residue.n_atoms ==  natoms:
            resID.append(res_id)
    #    print ('number of atoms in residue {}: {}'.format(res_id,residue.n_atoms))
    numberpolymers = len(resID)
    MoleculeResidueList = resID #range(0,numberpolymers) # important for Rg calculation
    print ("Unit cell:")
    print ("	{}".format(traj.unitcell_lengths[0])) # List of the unit cell on each frame
    print ('Number of frames:')
    print ("	{}".format(traj.n_frames))
    print ('Number of molecule types:')
    print ("	{}".format(traj.n_chains))
    print ('Number of molecules:')
    print ("	{}".format(traj.n_residues))
    print ('Number of atoms:')
    print ("	{}".format(traj.n_atoms))
    print ("Number of polymers:\n    {}".format(numberpolymers))
    print ("Number of atoms per chain:\n    {}".format(natoms))
    print ("Chain id of polymers:\n    {}".format(MoleculeResidueList))
    print ("Atom 1 coordinates:")
    print ('	{}'.format(traj.xyz[0][0]))
    # Get atom ends for Ree
    cnt = 0
    ReeAtomIndices = [] # remember atom ID - 1 is the atom index
    temp = []
    for i in range(numberpolymers*DOP):
        if cnt == 0:
            temp.append(i)
            cnt += 1
        elif cnt == (DOP-1):
            temp.append(i)
            ReeAtomIndices.append(temp)
            temp = []
            cnt = 0
        else:
            cnt += 1
            
    print(ReeAtomIndices)
    
    if Make_molecules_whole: # generates bonding list for each molecule
        bonds = []
        cnt = 0
        for i in range(numberpolymers*DOP):
            cnt += 1
            if cnt == (DOP):
                pass
                cnt = 0
            else:
                bonds_temp = [i,i+1]
                bonds.append(bonds_temp)
        if CnvDATA2PDB: # Automatically finds the bonds from the topology file
            bonds = None
        else:
            bonds = np.asarray(bonds,dtype=np.int32)
        traj.make_molecules_whole(inplace=True, sorted_bonds=bonds)
 
    if SaveCalculationsToFile == True: os.chdir(SaveFilename)
   
    ReeTimeseries = []
    RgTimeseries = []
    Ree_averages = []
    Rg_avg_stats = []
    Rg_averages = []
    Ree_avg_stats = []
    ''' Calculate Radius-of-gyration '''
    if CalcRg:
        print('Calculating radius-of-gyration...')
        # Compute the radius-of-gyration
        ElementDictionary ={
                    "carbon": 12.01,
                    "hydrogen": 1.008,
                    "oxygen": 16.00,
                    "nitrogen": 14.001,
                    "virtual site": 1.0,
                    "virtual_site": 1.0,
                    "sodium": "na+",
                    }

        Rg_list = []
        Rg_Ave_list = []

        for i,molecule in enumerate(MoleculeResidueList):
            atom_indices = traj.topology.select('resid {}'.format(i)) #and (resname UNL) or (resneme LEF) or (resname RIG)
            mass_list = []
            for index in atom_indices:
                temp = ElementDictionary[str(traj.topology.atom(index).element)]
                mass_list.append(temp)
            
            print ('Number of atoms in molecule {}'.format(i))
            print ('	{}'.format(len(atom_indices)))
            Rg = md.compute_rg(traj.atom_slice(atom_indices),np.asarray(mass_list))
            RgTimeseries.append(Rg)
            
            np.savetxt('Rg_out_mdtraj_molecule_{}.dat'.format((i)), Rg)

            stats_out = pymbar_statistics(Rg) #get time-averaged Rg of this molecule
            
            RgAvg = stats_out[0]
            RgVariance = stats_out[2]**2
            CorrTime = stats_out[5]
            Rg_averages.append([RgAvg,RgVariance,CorrTime])
            Rg_avg_stats.append([stats_out[6],stats_out[7],stats_out[8]])
            
            print ('The radius of gyration for molecule {} is:'.format(i))
            print ('	{0:2.4f} +/- {1:2.5f}'.format(RgAvg,np.sqrt(RgVariance)))
            
            ''' Plot the radius of gyration '''
            plt.plot(Rg, "k-")
            plt.xlabel('timestep')
            plt.ylabel('Radius-of-gryation')
            plt.savefig("Rg_molecule_{}.pdf".format(i),bbox_inches='tight')
            plt.close()

    if CalcRee:
        print('Calculating end-to-end distance...')
        for i,temp_pair in enumerate(ReeAtomIndices):
            EndEndDist = md.compute_distances(traj,atom_pairs=[temp_pair], periodic=False, opt=True)
            ReeTimeseries.append(EndEndDist)
            
            stats_out = pymbar_statistics(EndEndDist)
            
            ReeAvg = stats_out[0]
            ReeVariance = stats_out[2]**2
            CorrTime = stats_out[5]
            Ree_averages.append([ReeAvg,ReeVariance,CorrTime])
            Ree_avg_stats.append([stats_out[6],stats_out[7],stats_out[8]])
            
            print ('The End-end distance for molecule {} is:'.format(i))
            print ('	{0:2.4f} +/- {1:2.5f}'.format(ReeAvg,np.sqrt(ReeVariance)))
            
            ''' Plot the Ree '''
            plt.plot(EndEndDist, "k-")
            plt.xlabel('timestep')
            plt.ylabel('Ree')
            plt.savefig("Ree_{}.pdf".format(i),bbox_inches='tight')
            plt.close()

            np.savetxt('Ree_{}.dat'.format(i), EndEndDist)
    
    # Move backup to working directory
    if SaveCalculationsToFile == True: os.chdir("..")

    # Continue saving histograms to directories
    if SaveCalculationsToFile == True: os.chdir(SaveFilename)
    
    if CalcRee:
        ''' Histogram Ree '''
        Ree_temp = []
        for i in ReeTimeseries:
            Ree_temp.extend(i)
        Ree_data = np.asarray(Ree_temp)
#        HistogramRee(Ree_data, number_bins=25, DoBootStrapping=True, ShowFigures=False, NormHistByMax=True, 
#                            TrimRee=False, ReeCutoff=1.5, ReeMinimumHistBin=0., scale=1., gaussian_filter=False, sigma=2 )

        ''' Plot all the Ree '''
        for EndEndDist in ReeTimeseries:
            plt.plot(EndEndDist)
            plt.xlabel('timestep A.U.')
            plt.ylabel('Ree')
        plt.savefig("Ree_total.pdf",bbox_inches='tight')
        plt.close()

    if CalcRg:
        ''' Plot all the Rg '''
        for RgDist in RgTimeseries:
            plt.plot(RgDist)
            plt.xlabel('timestep A.U.')
            plt.ylabel('Rg')
        plt.savefig("Rg_total.pdf",bbox_inches='tight')
        plt.close()

    print ('****** TOTALS *******')

    ''' Calculate the ensemble averages '''
    stats_out = open('stats_out.data','w')
    if CalcRee:
        ReeTotal = 0
        Stats_ReeTotal = 0 
        ReeVarianceTotal = 0
        Stats_ReeVarTotal = 0 
        CorrTime = []
        Stats_CorrTime = []
        for index, Ree in enumerate(Ree_averages):
            ReeTotal = ReeTotal + Ree[0]
            ReeVarianceTotal = ReeVarianceTotal + Ree[1]
            CorrTime.append(Ree[2])
            
            #from stats.py script
            Stats_ReeTotal = Stats_ReeTotal + Ree_avg_stats[index][0]
            Stats_ReeVarTotal = Stats_ReeVarTotal + (Ree_avg_stats[index][1])**2
            Stats_CorrTime.append(Ree_avg_stats[index][2])
            

        ReeAverage = ReeTotal/len(Ree_averages)
        ReeStdErr  = np.sqrt(ReeVarianceTotal/len(Ree_averages))
        ReeAvgCorrTime = np.average(CorrTime)
        ReeCorrTimeStdErr = np.sqrt(np.var(CorrTime)/len(CorrTime))
        Stats_ReeAverage = Stats_ReeTotal/len(Ree_averages)
        Stats_StdErr = np.sqrt(Stats_ReeVarTotal/len(Ree_averages))
        Stats_AvgCorrTime = np.average(Stats_CorrTime)
        Stats_CorrTimeStdErr = np.sqrt(np.var(Stats_CorrTime)/len(CorrTime))
        print ('Total End-end distance average is: {0:4.4f} +/- {1:3.6f}'.format(ReeAverage,ReeStdErr))
        print ('Total End-end distance avg. correlation time: {0:5.4f} +/- {1:5.6f}'.format(ReeAvgCorrTime, ReeCorrTimeStdErr))
        print ('STATS: Total Ree distance avg is : {0:4.4f} +/- {1:3.6f}'.format(Stats_ReeAverage,Stats_StdErr))
        print ('STATS: Total Ree Corr. Time avg is : {0:4.4f} +/- {1:3.6f}'.format(Stats_AvgCorrTime,Stats_CorrTimeStdErr))
        stats_out.write('Total End-end distance average is: {0:4.4f} +/- {1:3.6f}\n'.format(ReeAverage,ReeStdErr))
        stats_out.write('Total End-end distance avg. correlation time: {0:5.4f} +/- {1:5.6f}\n'.format(ReeAvgCorrTime, ReeCorrTimeStdErr))
        stats_out.write('STATS: Total Ree distance avg is : {0:4.4f} +/- {1:3.6f}\n'.format(Stats_ReeAverage,Stats_StdErr))
        stats_out.write('STATS: Total Ree Corr. Time avg is : {0:4.4f} +/- {1:3.6f}\n'.format(Stats_AvgCorrTime,Stats_CorrTimeStdErr))
        
    if CalcRg:
        RgTotal = 0
        Stats_RgTotal = 0
        RgVarianceTotal = 0
        Stats_RgVarTotal = 0
        CorrTime = []
        Stats_RgCorrTime = []
        for index, Rg in enumerate(Rg_averages):
            RgTotal = RgTotal + Rg[0]
            RgVarianceTotal = RgVarianceTotal + Rg[1]
            CorrTime.append(Rg[2]) 
            
            #from stats.py script
            Stats_RgTotal = Stats_RgTotal + Rg_avg_stats[index][0]
            Stats_RgVarTotal = Stats_RgVarTotal + (Rg_avg_stats[index][1])**2
            Stats_RgCorrTime.append(Rg_avg_stats[index][2])
            
        RgAverage = RgTotal/len(Rg_averages)
        RgStdErr  = np.sqrt(RgVarianceTotal/len(Rg_averages))
        RgAvgCorrTime = np.average(CorrTime)
        RgCorrTimeStdErr = np.sqrt(np.var(CorrTime)/len(CorrTime))
        Stats_RgAverage = Stats_RgTotal/len(Rg_averages)
        Stats_RgStdErr = np.sqrt(Stats_RgVarTotal/len(Rg_averages))
        Stats_AvgRgCorrTime = np.average(Stats_RgCorrTime)
        Stats_RgCorrTimeStdErr = np.sqrt(np.var(Stats_RgCorrTime)/len(CorrTime))
        print ('Total Rg average is: {0:2.3f} +/- {1:2.5f}'.format(RgAverage,RgStdErr))
        print ('Total Rg avg. correlation time: {0:5.4f} +/- {1:5.6f}'.format(RgAvgCorrTime, RgCorrTimeStdErr))
        print ('STATS: Total Rg distance avg is : {0:4.4f} +/- {1:3.6f}'.format(Stats_RgAverage,Stats_RgStdErr))
        print ('STATS: Total Rg Corr. Time avg is : {0:4.4f} +/- {1:3.6f}'.format(Stats_AvgRgCorrTime,Stats_RgCorrTimeStdErr))
        stats_out.write('Total Rg average is: {0:2.3f} +/- {1:2.5f}\n'.format(RgAverage,RgStdErr))
        stats_out.write('Total Rg avg. correlation time: {0:5.4f} +/- {1:5.6f}\n'.format(RgAvgCorrTime, RgCorrTimeStdErr))
        stats_out.write('STATS: Total Rg distance avg is : {0:4.4f} +/- {1:3.6f}\n'.format(Stats_RgAverage,Stats_RgStdErr))
        stats_out.write('STATS: Total Rg Corr. Time avg is : {0:4.4f} +/- {1:3.6f}\n'.format(Stats_AvgRgCorrTime,Stats_RgCorrTimeStdErr))
        
        
    # Calculate alpha value
    if RgAverage != 0 and ReeAverage != 0:
        alpha = ReeAverage/RgAverage
        alpha_std = np.sqrt((1/RgAverage)**2*RgStdErr**2 + (ReeAverage/RgAverage**2)**2*ReeStdErr**2)
        Stats_alpha = Stats_ReeAverage/Stats_RgAverage
        Stats_alpha_std = np.sqrt((1/Stats_RgAverage)**2*Stats_RgStdErr**2 + (Stats_ReeAverage/Stats_RgAverage**2)**2*Stats_StdErr**2)
        print ('The alpha value: Ree/Rg is: {0:4.4f} +/- {1:4.4f}'.format(alpha,alpha_std))
        print ('STATS: The alpha value: Ree/Rg is: {0:4.4f} +/- {1:4.4f}'.format(Stats_alpha,Stats_alpha_std))
        stats_out.write('The alpha value: Ree/Rg is: {0:4.4f} +/- {1:4.4f}\n'.format(alpha,alpha_std))
        stats_out.write('STATS: The alpha value: Ree/Rg is: {0:4.4f} +/- {1:4.4f}\n'.format(Stats_alpha,Stats_alpha_std))
        
    # Move backup to working directory
    if SaveCalculationsToFile == True: os.chdir("..")

    with open(SysName+'_Statistics.dat','w') as g:
        g.write('#  Avg.    Pseudo-Var.    StdErr.     Corr.   Var.    StdErr.\n')
        if CalcRg:
            g.write('Rg        {0:8.4f}      {1:8.6f}      {2:8.6f}      {3:8.4f}      {4:8.6f}      {5:8.6f}\n'.format(Stats_RgAverage,Stats_RgVarTotal,Stats_RgStdErr,Stats_AvgRgCorrTime,np.var(Stats_RgCorrTime),Stats_RgCorrTimeStdErr))
        if CalcRee:
            g.write('Ree       {0:8.4f}      {1:8.6f}      {2:8.6f}      {3:8.4f}      {4:8.6f}      {5:8.6f}\n'.format(Stats_ReeAverage,Stats_ReeVarTotal,Stats_StdErr,Stats_AvgCorrTime,np.var(Stats_CorrTime),Stats_CorrTimeStdErr))

