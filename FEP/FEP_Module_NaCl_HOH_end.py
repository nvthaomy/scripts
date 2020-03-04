from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from parmed import gromacs
from parmed import unit as u
#from dcdreporter import DCDReporter
import parmed as pmd
import time
from openmmtools import *
import mdtraj as md
import mdtraj.reporters
import numpy as np
import stats_openmm as stats
from scipy.optimize import least_squares
import sys
"""Assumptions:
   If more than 1 residues, delete them with same frequency
   Topology sequence of system0 is n0 Na+, n0 Cl-, nw HOH
   Topology sequence of system1 is 
	1) inserting Na+ Cl-: n0 Na+, n0 Cl-, nw HOH, Na+, Cl-
        2) inserting HOH:     n0 Na+, n0 Cl-, nw HOH, HOH
   M.N: Modified for insertion in NPT (record box lengths in x y z and in every frames)
"""

''' ProgressBar Class '''

ProgressBarOn = sys.stdout.isatty()

class ProgressBar(object):
    def __init__(self, Text, Steps = 1, BarLen = 20, UpdateFreq = 1.):
        """Initializes a generic progress bar."""
        self.Text = Text
        self.Steps = Steps
        self.BarLen = BarLen
        self.UpdateFreq = UpdateFreq
        self.__LastTime = 0.
        self.__LastTime = time.time()
        self.__LastLen = 0
        self.Update(0)

    def Update(self, Step):
        """Updates the progress bar."""
        if time.time() - self.__LastTime > self.UpdateFreq:
            if not ProgressBarOn:
                return
            self.__LastTime = time.time()
            if self.BarLen == 0:
                s = "%s [%d]" % (self.Text, Step)
            else:
                Frac = float(Step) / (self.Steps + 1.e-300)
                n = int(self.BarLen * Frac + 0.5)
                n = max(min(n, self.BarLen), 0)
                s = "%s [" % self.Text
                s += "="*n + (self.BarLen-n)*" "
                s += "] %.2f%%" % (100.*Frac)
            self.__LastLen = len(s)
            s += "\r"
            sys.stdout.write(s)
            sys.stdout.flush()

    def Clear(self):
        """Clears text on this line."""
        if not ProgressBarOn:
            return
        sys.stdout.write(" "*self.__LastLen + "\r")
        sys.stdout.flush()

''' The FEP Class '''

class FEP:
    # constructor
    def __init__(self,method_,derivative_,traj_list_,states_list_,thermo_files_list_):
        ''' 
        
            method: 
                (1) Widom Insertions: WI (only for particle derivative)
                (2) Zwanzig: Z
                (3) Simple Overlap: SO
                (4) Non-optimal Bennett's: NOB
                (5) Optimal Bennett's: OB - recommended
            
            derivative: 
                (1) particle -> chemical potential 
                (2) volume -> pressure
                (3) area -> surface tension
                
            number states:
                (1) end-points: EP
                (2) single-state: SS - approximates end-point
                
            traj_list: 
                List of all mdtraj objects being used during reweighting. Assumed to 
                be in order 0,1, etc. (currenly only up to 2 should be there).
            
            states_list:
                OpenMM simulation objects being used during reweighting. Assumed to 
                be in order N,N+1,N-1 (currenly only up to 3 should be there). 
                
                if an end-points calc order is: 
                    (1) N,N+1
                    (2) V,V+dV
                
                if a single-state calc order is: 
                    (1) N,N+1,N-1
                    (2) V,V+dV,V-dV
                
            thermo_files_list:
                A list of list consisting of the potential energies for the reference states. 
        
        '''
        
        self.method = str(method_).lower()
        self.derivative = str(derivative_).lower()
        self.traj_list = traj_list_
        self.states_list = states_list_
        self.thermo_files_list = thermo_files_list_
        self.number_states = int(1) # default set to 1
        self.ConformerPositions = [] # list of possible conformers generated for insertions
        self.dU_0to1 = []
        self.dU_1to0 = []
        self.NumberDrawsPerFrame = int(1)
        self.NumberInsertionsPerFrame = int(1)
        self.NumberDeletionsPerFrame = int(1)
        self.ReRunRefState = False
        self.ThermalEnergy = float(2.479) # Thermal energy @ 298K Units KJ/mole
        self.Temperature = float(298.) # Kelvin
        self.BennettsConstant = float(0.) # Initial Guess for Bennett's Constant
        self.dF = float(0.) # the change in free energy (KJ/mole)
        self.dF_stderr = float(0.) # the standard error of change in free energy (KJ/mole)
        self.Perturbation = float(0.005) # the perturbation to apply to the system. For chemical potential, default is 1.
        self.IdealGasDOFs = int(1.) # for specifying the Degrees-of-freedom for the ideal gas calculation
        self.Perturbation_State0 = float(0.) # the valud of the perturbation parameter in state 0, the ref state
        self.Perturbation_State1 = float(0.) # the value of the perturbation parameter in state 1 (i.e. N+1, or V+dV)
        self.Perturbation_State2 = float(0.) # the valud of the perturbation parameter in state 2 (i.e. N-1, or V-dV)
        self.Pex = float(0.) # the excess pressure units atm
        self.Pex_stderr = float(0.) # the excess pressure standard error in units atm
        self.Pid = float(0.) # the ideal pressure units atm
        self.ResidueName = [] # Sets the residue names to be inserted/deleted, uses MDTraj selections
        self.ResidueIndiceList_Insertions = [] # the residue indices that can be chosen from traj for insertions
        self.ResidueIndiceList_Deletions = [] # the residue indices that can be deleted 
        self.ResidueNamesAndNumber_List = [] # the names and numbers of each residue in each traj in traj_list
        self.SurfaceTension = float(0.)
        self.SurfaceTension_stderr = float(0.)
        
        # Checks 
        
        if self.method == 'WI' and self.derivative == 'volume':
            print('ERROR: Widom Insertions only for derivative == particle')
            

    # methods
    def SetMethod(self,method_):
        ''' Set the method to run: OB, NOB, etc.'''
        self.method = str(method_).lower()
        
    def SetDerivative(self,derivative_):
        ''' Set the derivative to calculate: particle, volume, area, etc.'''
        self.derivative = str(derivative_).lower()
    
    def SetNumberOfDrawsPerFrame(self,NDrawsPF_):
        ''' Set the number of samples to draw from per frame 
            to generate molecule conformation library.'''
        self.NumberDrawsPerFrame = int(NDrawsPF_)
    
    def SetNumberInsertionsPerFrame(self,NIPF_):
        ''' Set the number of Insertions Per Frame '''
        self.NumberInsertionsPerFrame = int(NIPF_)
        
    def SetNumberDeletionsPerFrame(self,NDPF_):
        ''' Set the number of Deletions Per Frame '''
        self.NumberDeletionsPerFrame = (NDPF_)
        
    def SetReRunRefState(self,ReRunRefState_):
        ''' Set whether to rerun ref state, i.e. incase want to increase ewald tolerance during rerun '''
        self.ReRunRefState = ReRunRefState_
        
    def SetTemperature(self,Temperature_):
        ''' Set the temperature in units kelvin '''
        self.Temperature = float(Temperature_)
        self.ThermalEnergy = float(Temperature_*8.3145/1000.) # KJ/mole
        
    def SetBennettsConstant(self,BennettsConstant_):
        ''' Set the number of Deletions Per Frame '''
        self.BennettsConstant = float(BennettsConstant_)
        
    def SetPerturbation(self,Perturbation_):
        ''' Sets the Perturbation to apply '''
        self.Perturbation = float(Perturbation_)
        
    def SetIdealGasDOFs(self,IdealGasDOFs_):
        ''' Sets the degrees of freedom for the corresponding ideal gas'''
        self.IdealGasDOFs = int(IdealGasDOFs_)
        
    def SetResidueName(self,ResidueName_):
        ''' Sets the residue name to be used in insertions and deletions,
            from the topology file loaded in. Used MDTraj selection commands to do this. '''
        self.ResidueName = ResidueName_
        
    def SetResidueIndiceList_Insertions(self,ResidueIndiceList_):
        ''' A list of the residues in the topology file that can be chosen for Insertions. '''
        self.ResidueIndiceList_Insertions = ResidueIndiceList_
        
    def SetResidueIndiceList_Deletions(self,ResidueIndiceList_):
        ''' A list of the residues in the topology file that can be chosen for Insertions. '''
        self.ResidueIndiceList_Deletions = ResidueIndiceList_
        
    def GetResidueNames(self):
        ''' Returns a list with the names and number of the residues, and number atoms in each traj object.  '''
        ResidueNamesANDNumber_List = []
        for i,traj in enumerate(self.traj_list):
            temp = []
            temp_name = []
            temp_atoms = [] # number of atoms in a residue of type name
            residues = traj.topology.residues
            for residue in residues:
                if residue.name in temp_name: # do not add if already in list
                    pass
                else: # add new residue type to list
                    temp_name.append(str(residue.name))
                    temp_atoms.append(int(residue.n_atoms))
            
            for j,res_name in enumerate(temp_name): # find the number of residues
                number_resname = len(traj.topology.select("resname '{}'".format(str(res_name))))
                temp.extend([str(res_name),int(int(number_resname)/int(temp_atoms[j])),int(temp_atoms[j])])
                
            ResidueNamesANDNumber_List.append(temp)
        
        self.ResidueNamesAndNumber_List = ResidueNamesANDNumber_List
        
        return ResidueNamesANDNumber_List
        
    def CalcPex(self):
        ''' Calculates the excess pressure, units atm '''
        dF = self.dF # units KJ/mole
        
        dF_stderr = self.dF_stderr # units KJ/mole
        dV = float(self.Perturbation_State1-self.Perturbation_State0) # dV in units nm^3
        
        Nav = 6.022E23 # atoms/mole
        Pa2Atm = 101325. # number of Pascals in 1 atm
        
        Pex = -1*dF*1000./Nav/(1E-9)**3/dV/Pa2Atm
        Pex_stderr = (1000./Nav/(1E-9)**3/dV/Pa2Atm)*dF_stderr
        
        self.Pex = Pex
        self.Pex_stderr = Pex_stderr
        
    def CalcPid(self):
        ''' Calculates the ideal pressure, units atm '''
        NDOF= float(self.IdealGasDOFs) 
        V_State0 = float(self.Perturbation_State0) # in units nm^3
        Nav = 6.022E23 # atoms/mole
        Pa2Atm = 101325. # number of Pascals in 1 atm
        
        Pid = self.ThermalEnergy*NDOF*1000./Nav/(1E-9)**3/V_State0/Pa2Atm
        
        self.Pid = Pid    
        
    def CalcSurfaceTension(self):
        ''' Calculates the surface tension '''
        dF = self.dF # units KJ/mole
        dF_stderr = self.dF_stderr # units KJ/mole
        dA = float(self.Perturbation_State1-self.Perturbation_State0) # dA in units nm^2
        
        Nav = 6.022E23 # atoms/mole
        
        SurfTen = dF*1000./Nav/(1E-9)**2/2./dA # divide by factor of 2 to acocunt for both interfaces
        SurfTen_stderr = (1000./Nav/(1E-9)**2/dA/2.)*dF_stderr
        
        self.SurfaceTension = SurfTen
        self.SurfaceTension_stderr = SurfTen_stderr
        
    def BennettsWeightingFunction(self,x,beta,C):
        
        weights = 1./np.cosh(beta*(x-C)/2.)
        
        return weights
        
    def BennettsDeltaFreeEnergy(self,BennettsConstant,dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta):
        
        C = BennettsConstant
        
        weights_0to1 = self.BennettsWeightingFunction(dU_0to1,beta,C)
        weights_1to0 = self.BennettsWeightingFunction((-1*dU_1to0),beta,C)
        
        temp = np.multiply(weights_0to1,ExpdU_0to1)
        g = np.average(temp)
        g_stddev = np.std(temp)
        g_count  = temp.shape[0]
        g_stderr = g_stddev/np.sqrt(g_count)
        
        temp = np.multiply(weights_1to0,ExpdU_1to0)
        h = np.average(temp)
        h_stddev = np.std(temp)
        h_count  = temp.shape[0]
        h_stderr = h_stddev/np.sqrt(h_count)
        
        dF = -1./beta*np.log(g/h)
        
        stderr = np.sqrt(((1./beta/g)**2*g_stderr**2+(1./beta/h)**2*h_stderr**2))
        
        return dF,stderr
        
    def DeltaFreeEnergy(self,ExpdU_0to1,ExpdU_1to0,beta):
        
        temp = ExpdU_0to1
        g = np.average(temp)
        g_stddev = np.std(temp)
        g_count  = temp.shape[0]
        g_stderr = g_stddev/np.sqrt(g_count)
        
        temp = ExpdU_1to0
        h = np.average(temp)
        h_stddev = np.std(temp)
        h_count  = temp.shape[0]
        h_stderr = h_stddev/np.sqrt(h_count)
        
        dF = -1./beta*np.log(g/h)
        
        stderr = np.sqrt(((1./beta/g)**2*g_stderr**2+(1./beta/h)**2*h_stderr**2))
        
        return dF,stderr
        
    def BennettsCostFunction(self,x,dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta):
        
        BennettsConstant = x
        
        dF,stderr = self.BennettsDeltaFreeEnergy(BennettsConstant,dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta)
        residual = (dF - BennettsConstant) # These ought to be equal for minimum variance
        
        return residual
        
    def Widom_Insertions():
        ''' Uses Widom Insertions to estimate the chemical potential. '''
        
    def Zwangzig():
        ''' The Zwangzig relation is essential Widom's Insertion method, but for any perturbation. '''
        
    def Simple_Overlap():
        ''' Simple, two state overlap estimate, as compared to Bennett's method, the weighting function here is set to 1.'''
        
        dU_0to1 = np.asarray(self.dU_0to1)[:,2] # in KJ/mole
        dU_1to0 = np.asarray(self.dU_1to0)[:,2] # in KJ/mole
        beta = 1./self.ThermalEnergy # inverse KJ/mole
        
        ExpdU_0to1 = np.exp(-1*beta*dU_0to1/2.)
        ExpdU_1to0 = np.exp(-1*beta*dU_1to0/2.)
        
        dF,dF_stderr = self.DeltaFreeEnergy(ExpdU_0to1,ExpdU_1to0,beta)
        
        # set the change in free energy and the standard error
        self.dF = dF
        self.dF_stderr = dF_stderr        
        
    def NonOptimal_Bennetts():
        ''' Non-optimal Bennett's Method, i.e. constant==specified by user '''
            
        dU_0to1 = np.asarray(self.dU_0to1)[:,2] # in KJ/mole
        dU_1to0 = np.asarray(self.dU_1to0)[:,2] # in KJ/mole
        beta = 1./self.ThermalEnergy # inverse KJ/mole
        
        ExpdU_0to1 = np.exp(-1*beta*dU_0to1/2.)
        ExpdU_1to0 = np.exp(-1*beta*dU_1to0/2.)
        
        OptimalBennettsConstant = 0.
        
        self.BennettsConstant = float(OptimalBennettsConstant)
        
        dF,dF_stderr = self.BennettsDeltaFreeEnergy(OptimalBennettsConstant,dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta)
        
        # set the change in free energy and the standard error
        self.dF = dF
        self.dF_stderr = dF_stderr
        
    def Optimal_Bennetts(self):
        ''' Performs FEP using Bennett's Method, the best method. '''
        C_init = self.BennettsConstant
            
        dU_0to1 = np.asarray(self.dU_0to1)[:,2] # in KJ/mole
        dU_1to0 = np.asarray(self.dU_1to0)[:,2] # in KJ/mole
        beta = 1./self.ThermalEnergy # inverse KJ/mole
        
        ExpdU_0to1 = np.exp(-1*beta*dU_0to1/2.)
        ExpdU_1to0 = np.exp(-1*beta*dU_1to0/2.)
        
        opt = least_squares(self.BennettsCostFunction,self.BennettsConstant,args=(dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta))

        OptimalBennettsConstant = opt.x
        
        self.BennettsConstant = float(OptimalBennettsConstant)
        
        dF,dF_stderr = self.BennettsDeltaFreeEnergy(OptimalBennettsConstant,dU_0to1,dU_1to0,ExpdU_0to1,ExpdU_1to0,beta)
        
        # set the change in free energy and the standard error
        self.dF = dF
        self.dF_stderr = dF_stderr
    
    def ReduceCoordinates(self,xyz,boxvectors):
        ''' Puts all the coordinates in reduced units, i.e. x_r = x/L_x, etc. '''
        Lbox = [boxvectors[0][0],boxvectors[1][1],boxvectors[2][2]]
        xyz_red = np.divide(xyz,Lbox,dtype='float32')
        
        return xyz_red
        
    def UnReduceCoordinates(self,xyz_red,new_boxvectors):
        ''' rescales all reduced coordinates in a frame to new box coordinates '''
        Lbox = [new_boxvectors[0][0],new_boxvectors[1][1],new_boxvectors[2][2]]
        xyz = np.multiply(xyz_red,Lbox,dtype='float32')
        
        return xyz
    
    def GenerateMoleculeConfigsInsertions(self,traj,DrawsPerFrame,NumberInsertionsPerFrame):
        '''
            Generates Molecule Configurations for insertions from trajectory frames.
            
            traj: mdtraj object (loaded not iterator)
            DrawsPerFrame: number of random molecule selections per frame of traj
        '''

        NumberFrames             = traj.n_frames
        NumberMolecules         = traj.n_residues # assumed to be the number of molecules
        DrawPerFrame             = DrawsPerFrame
        BoxVecs = traj.unitcell_vectors #box vectors of all frames
        BoxLs = [Vec.diagonal() for Vec in BoxVecs] #box length in x y z of all frames
       
        temp_top = traj.topology.copy()
        
        temp_resid_lists = []
        if self.ResidueName == []: # no residue set
            res_name = temp_top.residue(0).name # The model residue that is to be inserted
            for resid,residue in enumerate(temp_top.residues):
                if res_name == residue.name:
                    temp_resid_list.append(resid)
            model_res = temp_resid_list # indices of residues of type res_name
            NumberMolecules = len(model_res)
        else: 
            res_names = self.ResidueName # the model residue is set
            for res_name in res_names:
                temp_resid_list = []
                for resid,residue in enumerate(temp_top.residues):
                    if res_name in residue.name:
                        temp_resid_list.append(resid)
                
                temp_resid_lists.append(temp_resid_list)
            model_res = temp_resid_lists #list of lists of indices of residues of type res_names
            NumberMolecules = len(model_res[0])
            
#        self.SetResidueIndiceList_Insertions([x in y for y in model_res]) # for record keeping
        
        if NumberMolecules < 10:
            print('WARNING: There are only {} molecules to select from for insertions. Is this correct?'.format(int(NumberMolecules)))
        
        #atoms_in_residue = len(model_res)
        atoms_masses_list = [] #list of lists of atom masess in residues to be inserted
        for i, res_name in enumerate(res_names):      
            atoms_masses = []
            for atom in temp_top.residue(model_res[i][0]).atoms: # add atoms to temp_residue
                t_name         = atom.name
                t_element     = atom.element
                t_mass         = atom.element.mass
                atoms_masses.append(t_mass)
            atoms_masses_list.append(atoms_masses)
             
        if DrawPerFrame == 0:
            DrawPerFrame = 1
        elif DrawPerFrame > NumberMolecules:
            DrawPerFrame = NumberMolecules
            print('WARNING: DrawPerFrame > NumberMolecules ; setting to NumberMolecules.'.format(NumberMolecules))
        
        # Generate random integers for selecting molecules 
        MoleculeIndices2DrawList = []
        for resid_list in model_res:
            MoleculeIndices2Draw    =  np.random.randint(0,len(resid_list),size=(DrawPerFrame*NumberFrames))
            MoleculeIndices2Draw    = [resid_list[i] for i in MoleculeIndices2Draw]
            MoleculeIndices2DrawList.append(MoleculeIndices2Draw) 
         
        LibraryPositions_Lists    = []
        LibraryCOM_Lists            = []
        masses_list                     = atoms_masses_list
        
        PBar = ProgressBar('Drawing Conformers:', Steps = (NumberFrames-1), BarLen = 20, UpdateFreq = 1.)
        
        for j, MoleculeIndices2Draw in enumerate(MoleculeIndices2DrawList):
            masses = masses_list[j]
            LibraryPositions_List    = []
            LibraryCOM_List          = []
            for ind, frame in enumerate(traj): # iterate through trajectory and get molecule positions
                PBar.Update(ind)
                for i in range(DrawPerFrame):
                    # select atom indices belonging to molecule (i.e. to the residue)
                    temp_atom = frame.topology.select('resid {}'.format(MoleculeIndices2Draw[(i+DrawPerFrame*ind)]))
                    temp_pos = frame.xyz[0][temp_atom]
                    temp_com = np.multiply(np.transpose(temp_pos),masses)
                    temp_sum_com = np.sum(temp_com,axis=1)
                    temp_sum_mass = np.sum(masses)
                    temp_com = np.divide(temp_sum_com,temp_sum_mass)
                    temp_pos = np.subtract(temp_pos,temp_com) # Subtract off the COM, we will randomly regenerate this
                    LibraryPositions_List.append(temp_pos)
                    # COM currently not used again outside this function, but still create a library
                    LibraryCOM_List.append(temp_com)
            LibraryPositions_Lists.append(LibraryPositions_List)
            LibraryCOM_Lists.append(LibraryCOM_List)        
        PBar.Clear()
        
        #making BoxLs matrix to match NewCOM
        BoxLs_N = []
        for i in range(NumberFrames):
            Ls = BoxLs[i]
            BoxLs_N.append([Ls]*NumberInsertionsPerFrame)
        BoxLs = np.array(BoxLs_N)
        if BoxLs.shape != (NumberFrames,NumberInsertionsPerFrame,3):
            Exception('Size of box length matrix BoxLs does not match (NumberFrames,NumberInsertionsPerFrame,3)')

        NewCOMs = []
        NewPositions_Lists = []
        for i, LibraryPositions_List in enumerate(LibraryPositions_Lists):
            NewPositions_List = []
            #Randomly Generate new COM's
            NewCOM = BoxLs * np.random.uniform(size=(NumberFrames,NumberInsertionsPerFrame,3))
            NumberConformationsInLibrary = len(LibraryPositions_List)
            
            NewCOMs.append(NewCOM)
        
            #Randomly select new samples from library
            LibrarySelection = np.random.randint(0,(NumberConformationsInLibrary),size=(NumberFrames,NumberInsertionsPerFrame))
    
            # Generate New Positions
            for i in range(NumberFrames):
                temp = []
                for j in range(NumberInsertionsPerFrame):
                    temp1 = np.add(LibraryPositions_List[LibrarySelection[i][j]],NewCOM[i][j])
                    temp.append(temp1)
                NewPositions_List.append(temp)
            NewPositions_Lists.append(NewPositions_List)
         
        #to check that NewPositions_Lists has correct dimension
        try:
            test = np.array(NewPositions_Lists)
        except:
            Exception('Lengths of position lists are not the same between residues')

        return NewPositions_Lists
        
    def PerformInsertions(self,traj,simulation0,simulation1,NewPositions_Lists,PotEne_State0_List,ReRunRefState):
        ''' -Performs the Insertions to Calculate dU for state 0 to 1.
            -RefSate is always assumed to be state0.
        
        '''
        NumberFrames = traj.n_frames
        atom_indices = range(traj.n_atoms) 
        nAtoms0 = traj.n_atoms
        PotEne_State0 = [] # used if rerunning the ref. state as well.
        PotEne_Data = [] 
        waterId = None        
        CLId = None 
        PBar = ProgressBar('Insertion Progress:', Steps = (NumberFrames-1), BarLen = 20, UpdateFreq = 1.)
        
        # Add a residue for particle inserting
        temp_top0 = traj.topology.copy()
        temp_top = traj.topology.copy() 
        n_chains = temp_top.n_chains
        
        model_res_List = []
        if self.ResidueName == []: # no residue set
            model_res_List = [temp_top.residue(0)] # The model residue that is to be inserted
        else: 
            res_names = self.ResidueName # the model residue is set
            for res_name in res_names:
                for resid,residue in enumerate(temp_top.residues):
                    if residue.name == res_name:
                        model_res_List.append(temp_top.residue(resid))
                        break
        #get index of first water molecule
        for resid,residue in enumerate(temp_top.residues):
            if residue.name == 'HOH' and waterId == None:
                waterId = resid
                print('index of first water molecule is {}'.format(waterId))
            if residue.name == 'CL-' and CLId == None:
                CLId = resid
                print('index of first Cl-  molecule is {}'.format(CLId))
        if model_res_List == []: # check if blank
            print('ERROR: No model residue found.')
            
        atoms_added_List = []
        atoms_added_indices_List = []
        for i, model_res in enumerate(model_res_List):

            temp_chain = temp_top.add_chain() # Create a new chain
            temp_residue = temp_top.add_residue(model_res.name,temp_chain) # Create a new residue on temp_chain
            atoms_added = 0
            atoms_added_indices = []
            for atom in temp_top.residue(model_res.index).atoms: # add atoms to temp_residue
                t_name = atom.name
                t_element = atom.element    
                t_atom = temp_top.add_atom(t_name, t_element, temp_residue, serial=None)
                atoms_added += 1 
                atoms_added_indices.append(t_atom.index)
            atoms_added_List.append(atoms_added)
            atoms_added_indices_List.append(atoms_added_indices)
 
        NumAtoms = temp_top.n_atoms # Number of atoms after one insertion
        xyz = np.zeros((1,NumAtoms,3))
        tot_atoms_added = sum(atoms_added_List)

        for i in range(NumberFrames):
            PBar.Update(i)
            time_start = time.time() # time how long it takes
            # update the initial positions, positions for additional NaCl pair is still 0
            xyz[0][:nAtoms0] = traj[i].xyz[0]
            box =  traj[i].unitcell_vectors[0]
            
            if ReRunRefState:
                simulation0.context.setPeriodicBoxVectors(box[0],box[1],box[2])
                simulation0.context.setPositions(traj[i].xyz[0])
                simulation0.context.applyConstraints(1e-12)
                state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
                temp_PotEne_0 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            else: 
                temp_PotEne_0 = PotEne_State0_List[i]
#            print(temp_PotEne_0) 
            PotEne_State0.append(temp_PotEne_0)
            temp_traj0 = md.Trajectory(traj[i].xyz,temp_top0)
            temp_traj0.save_pdb('traj0_Insert.pdb') 

            ''' Get N+1 Potential Energy '''
            for j in range(self.NumberInsertionsPerFrame):
                for k, atoms_added in enumerate(atoms_added_List):
                    NewPositions_List = NewPositions_Lists[k]
                    if self.ResidueName[k] in ['Na+','NA+']:
                        t_index = nAtoms0 
                        xyz[0][t_index] = NewPositions_List[i][j][0]
                    elif self.ResidueName[k] in ['CL-','Cl-']:
                        t_index = nAtoms0+1
                        if nAtoms0+1 != NumAtoms: #Check if number of atoms is correct
                            Exception('Index of inserted atoms does not match total number of atoms')
                        xyz[0][t_index] = NewPositions_List[i][j][0]
                    elif self.ResidueName[k] in ['WAT','HOH']:
                        for ind in range(atoms_added):
                            t_index = int(atoms_added-ind)
                            xyz[0][(-1*t_index)] = NewPositions_List[i][j][ind]
                
                temp_traj = md.Trajectory(xyz,temp_top)
                
                #temp_traj = md.Trajectory(xyz,traj_System2.topology)
                temp_traj.save_pdb('traj1_Insert.pdb')
                
                # Get energy from state1 using new positions
                simulation1.context.setPeriodicBoxVectors(box[0],box[1],box[2])
                simulation1.context.setPositions(temp_traj.xyz[0])
                simulation1.context.applyConstraints(1e-12)
                state1 = simulation1.context.getState(getEnergy=True,enforcePeriodicBox=False)
                temp_PotEne_1 = state1.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                
                dU = temp_PotEne_1-temp_PotEne_0 # U1-U0
                

                PotEne_Data.append([temp_PotEne_0, temp_PotEne_1, dU])
                flog = open('ReRun_PotEne_State0to1_Insertions.log','a')
                flog.write('{},  {},  {},  {},  {}\n'.format(i,j,temp_PotEne_0, temp_PotEne_1, dU))
                flog.close()
            
            time_end = time.time()
            #print('Time: {0:3.5f}'.format(time_end-time_start))
        
        PBar.Clear()
        
        return PotEne_Data, PotEne_State0
    
    def PerformDeletions(self,traj,simulation0,simulation1,NumberDeletionAttemptsPerFrame,PotEne_State0_List,ReRunRefState):
        ''' -Performs Deletions to Calculate dU for state 0 to 1; i.e. N to N-1. IF an end-point simulation is being 
                conducted, with 2 reference simulations, the deletions occur for N+1 to N. 
            -RefSate is always assumed to be state0.
        '''
        
        NumberFrames = traj.n_frames
        BoxL = traj[0].unitcell_vectors[0][0][0]
        NumberMolecules = traj[0].n_residues
        PotEne_Data = []
        PBar = ProgressBar('Deletion Progress:', Steps = (NumberFrames-1), BarLen = 20, UpdateFreq = 1.)
        
        # Create list of indices of residues that can be deleted
        temp_top0 = traj.topology.copy()
        temp_top = traj.topology.copy()
        model_res_List = []     
        temp_resid_lists = []
        if self.ResidueName == []: # no residue set
            res_name = temp_top.residue(0).name # The model residue that is to be deleted
            for resid,residue in enumerate(temp_top.residues):
                if res_name == residue.name:
                    temp_resid_list.append(resid)
            model_res = temp_resid_list # indices of residues of type res_name
            NumberMolecules = len(model_res)
        else: 
            res_names = self.ResidueName # the model residue is set
            for i,res_name in enumerate(res_names):
                temp_resid_list = []
                for resid,residue in enumerate(temp_top.residues):
                    if res_name == residue.name:
                        temp_resid_list.append(resid)
                temp_resid_lists.append(temp_resid_list)
            model_res = temp_resid_lists # indices of residues of types res_names
            print('NumberMolecules {}'.format(NumberMolecules))
            NumberMolecules = len(model_res[0])
        
        self.SetResidueIndiceList_Deletions(model_res) # for record keeping
        
        if NumberMolecules< 10:
            print('WARNING: There are only {} molecules to select from for deletions. Is this correct?'.format(int(NumberMolecules)))

        atoms_masses_list=[]
        for i,resid in enumerate(model_res):
            atoms_in_residue = len(resid)
            atoms_masses = []            
            for atom in temp_top.residue(resid[0]).atoms: # add atoms to temp_residue
                t_name         = atom.name
                t_element     = atom.element
                t_mass         = atom.element.mass
                atoms_masses.append(t_mass)
            atoms_masses_list.append(atoms_masses)

        if NumberDeletionAttemptsPerFrame == 0:
            NumberDeletionAttemptsPerFrame = 1
        elif NumberDeletionAttemptsPerFrame > NumberMolecules:
            NumberDeletionAttemptsPerFrame = NumberMolecules
            print('WARNING: NumberDeletionAttemptsPerFrame > NumberMolecules ; setting to NumberMolecules {}.'.format(NumberMolecules))
        
        MoleculeIndices2DeleteList = []
        RandomSelect = False
        for i,resid in enumerate(model_res):
            # Select which residues to delete
            MoleculeIndices2Delete = []
            if RandomSelect: # Not Recommended
                temp = np.random.randint(0,(NumberMoleculesList[0]),size=(NumberFrames,NumberDeletionAttemptsPerFrame))
                for j in temp: 
                    temp_temp = [resid[i] for i in j]
            
                MoleculeIndices2Delete.append(temp_temp)
            
            else:
                temp = np.linspace(0,NumberDeletionAttemptsPerFrame-1,NumberDeletionAttemptsPerFrame)
                MoleculeIndices2Delete = [resid[int(i)] for i in temp]
            MoleculeIndices2DeleteList.append(MoleculeIndices2Delete)

        for i in range(NumberFrames):
            
            PBar.Update(i)
            time_start = time.time() # time how long it takes
            box =  traj[i].unitcell_vectors[0]
 
            temp_frame = traj[i].xyz[0]
            if ReRunRefState:
                simulation0.context.setPeriodicBoxVectors(box[0],box[1],box[2])
                simulation0.context.setPositions(traj[i].xyz[0])
                simulation0.context.applyConstraints(1e-12)
                state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
                temp_PotEne_0 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            else: 
                temp_PotEne_0 = PotEne_State0_List[i]            
            
            temp_traj0 = md.Trajectory(traj[i].xyz,temp_top0)
            temp_traj0.save_pdb('traj0_Delete.pdb')
            temp_top = traj[i].topology
#            print(temp_PotEne_0)    

            for j in range(NumberDeletionAttemptsPerFrame):
                select_expression = ''
                for k, MoleculeIndices2Delete in enumerate(MoleculeIndices2DeleteList):
                    if RandomSelect:
                        s = '(resid != {})'.format(MoleculeIndices2Delete[i][j]) 
                    else:
                        s = '(resid != {})'.format(int(MoleculeIndices2Delete[j]))

                    if len(select_expression) == 0:
                        select_expression += s
                    else:
                        select_expression += ' and {}'.format(s)
                residueSelect = temp_top.select(select_expression)
                temp_traj2 = traj[i].atom_slice(residueSelect)
                temp_traj2.save_pdb('traj2_Delete.pdb')

                xyz = traj[i].atom_slice(residueSelect).xyz # Select all other atoms, other than one deleting
                simulation1.context.setPeriodicBoxVectors(box[0],box[1],box[2])
                simulation1.context.setPositions(xyz[0])
                simulation1.context.applyConstraints(1e-12)
                state = simulation1.context.getState(getEnergy=True,enforcePeriodicBox=False)
                temp_PotEne_1 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                
                dU = (temp_PotEne_1-temp_PotEne_0) 
                
                PotEne_Data.append([temp_PotEne_0, temp_PotEne_1, dU])
                flog = open('ReRun_PotEne_State1to0_Deletions.log','a')
                flog.write('{},  {},  {},  {},  {}\n'.format(i,j,temp_PotEne_0, temp_PotEne_1, dU))
                flog.close()    
        
        PBar.Clear()
        
        return PotEne_Data
    
    def ApplyVolumePerturbation(self,traj,simulation0,PotEne_State0_List,ReRunRefState):
        ''' -Performs a volume perturbation 
            -RefSate is always assumed to be state0.
        '''
        NumberFrames = traj.n_frames
        BoxL = traj[0].unitcell_vectors[0][0][0]
        box =  traj[0].unitcell_vectors[0]        
        NumberMolecules = traj[0].n_residues
        boxvolume = simulation0.context.getState().getPeriodicBoxVolume()/nanometer**3 # units nanometers^3, volume of state 0
        periodicboxvectors = simulation0.context.getState().getPeriodicBoxVectors()/nanometer # get periodic box vectors, in nanometer
        
        V_expansion = boxvolume*(1.+float(self.Perturbation))
        V_contraction = boxvolume*(1.-float(self.Perturbation))
        
        Lbox_exp = V_expansion**(1./3.)
        Lbox_cont = V_contraction**(1./3.)
        
        PBV_expansion = [[Lbox_exp,0.,0.],[0.,Lbox_exp,0.],[0.,0.,Lbox_exp]]*nanometer
        PBV_contraction = [[Lbox_cont,0.,0.],[0.,Lbox_cont,0.],[0.,0.,Lbox_cont]]*nanometer
        
        self.Perturbation_State0 = float(boxvolume) # the valud of the perturbation parameter in state 0, the ref state
        self.Perturbation_State1 = float(V_expansion) # the value of the perturbation parameter in state 1 (i.e. N+1, or V+dV)
        self.Perturbation_State2 = float(V_contraction) # the valud of the perturbation parameter in state 2 (i.e. N-1, or V-dV)
        
        
        PotEne_Data_Expansion = []
        PotEne_Data_Contraction = []
        PotEne_State0 = [] # used if rerunning the ref. state as well. 
        
        PBar = ProgressBar('Volume Perturbation Progress:', Steps = (NumberFrames-1), BarLen = 20, UpdateFreq = 1.)
        
        for i in range(NumberFrames):
            
            PBar.Update(i)
            
            xyz = traj[i].xyz[0]
            xyz_red = self.ReduceCoordinates(xyz,periodicboxvectors) # get reduced coordinates
            
            time_start = time.time() # time how long it takes
            
            if ReRunRefState:
                simulation0.context.setPeriodicBoxVectors(box[0],box[1],box[2])
                simulation0.context.setPositions(traj[i].xyz[0])
                simulation0.context.applyConstraints(1e-12)
                state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
                temp_PotEne_0 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            else: 
                temp_PotEne_0 = PotEne_State0_List[i]            
            
            PotEne_State0.append(temp_PotEne_0)
            
            
            
            
            
            # DO EXPANSION
            xyz_exp = self.UnReduceCoordinates(xyz_red,PBV_expansion/nanometer)
            simulation0.context.setPeriodicBoxVectors(PBV_expansion[0],PBV_expansion[1],PBV_expansion[2])
            simulation0.context.setPositions(xyz_exp)
            simulation0.context.applyConstraints(1e-12)
            state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
            temp_PotEne_1 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            
            dU_0to1 = (temp_PotEne_1-temp_PotEne_0) 
            
            PotEne_Data_Expansion.append([temp_PotEne_0, temp_PotEne_1, dU_0to1])
            flog = open('ReRun_PotEne_State0to1_VolExpansion.log','a')
            flog.write('{},  {},  {},  {}\n'.format(i,temp_PotEne_0, temp_PotEne_1, dU_0to1))
            flog.close()    
        
            # DO CONTRACTION
            xyz_cont = self.UnReduceCoordinates(xyz_red,PBV_contraction/nanometer)
            temp_traj = md.Trajectory(xyz_cont,traj.topology)
            simulation0.context.setPeriodicBoxVectors(PBV_contraction[0],PBV_contraction[1],PBV_contraction[2])
            simulation0.context.setPositions(xyz_cont)
            simulation0.context.applyConstraints(1e-12)
            state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
            temp_PotEne_2 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            
            dU_1to0 = (temp_PotEne_2-temp_PotEne_0) 
            
            PotEne_Data_Contraction.append([temp_PotEne_0, temp_PotEne_2, dU_1to0])
            flog = open('ReRun_PotEne_State0to1_VolContraction.log','a')
            flog.write('{},  {},  {},  {}\n'.format(i,temp_PotEne_0, temp_PotEne_2, dU_1to0))
            flog.close()    
        
        PBar.Clear()
    
        return PotEne_Data_Expansion, PotEne_Data_Contraction, PotEne_State0
    
    def ApplyAreaPerturbation(self,traj,simulation0,PotEne_State0_List,ReRunRefState):
        ''' -Performs an area perturbation 
            -RefSate is always assumed to be state0.
        '''
        NumberFrames = traj.n_frames
        BoxL = traj[0].unitcell_vectors[0][0][0]
        box =  traj[0].unitcell_vectors[0]        
        NumberMolecules = traj[0].n_residues
        boxvolume = simulation0.context.getState().getPeriodicBoxVolume()/nanometer**3 # units nanometers^3, volume of state 0
        periodicboxvectors = simulation0.context.getState().getPeriodicBoxVectors()/nanometer # get periodic box vectors, in nanometer
        
        #Find Dimension with the largest length 
        Lmax = BoxL
        dim_max = 0
        for dim, vector in enumerate(box):
            if vector[dim]>Lmax:
                Lmax = float(vector[dim]) # units nanometers
                dim_max = dim # 0 = x; 1 = y ; 2 = z
                
        print('Dimension {} is the max with length {} nms.'.format(dim_max,Lmax))
    
        #scaling1 = expanding area, scaling2 = shrinking area
        dAfrac = self.Perturbation
        scaleZ1     = 1.0/(1.0 + dAfrac)
        scaleXY1    = (1.0 + dAfrac)**0.5
        scaleZ2     = 1.0/(1.0 - dAfrac)
        scaleXY2    = (1.0 - dAfrac)**0.5
        
        
        Area = float(boxvolume/Lmax) # units nm**2
        
        A_expansion = Area*(1.+float(self.Perturbation))
        A_contraction = Area*(1.-float(self.Perturbation))
        
        # maintain constant volume
        Lmax_expansion = boxvolume/(A_expansion)
        Lmax_contractrion = boxvolume/(A_contraction)
        
        Lbox_exp = A_expansion**(1./2.) # units nm
        Lbox_cont = A_contraction**(1./2.) # units nm
        
        
        PBA_expansion = []
        PBA_contraction = []
        for dim in range(3): # create new box dimensions
            temp_exp = [0.,0.,0.]
            temp_cont = [0.,0.,0.]
            if dim == dim_max:
                temp_exp[dim] = Lmax_expansion
                temp_cont[dim] = Lmax_contractrion
            else:
                temp_exp[dim] = Lbox_exp
                temp_cont[dim] = Lbox_cont
            PBA_expansion.append(temp_exp)
            PBA_contraction.append(temp_cont)
                
        PBA_expansion = PBA_expansion*nanometer
        PBA_contraction = PBA_contraction*nanometer
        
        self.Perturbation_State0 = float(Area) # the valud of the perturbation parameter in state 0, the ref state
        self.Perturbation_State1 = float(A_expansion) # the value of the perturbation parameter in state 1 (i.e. N+1, or V+dV)
        self.Perturbation_State2 = float(A_contraction) # the valud of the perturbation parameter in state 2 (i.e. N-1, or V-dV)
        
        
        PotEne_Data_Expansion = []
        PotEne_Data_Contraction = []
        PotEne_State0 = [] # used if rerunning the ref. state as well. 
        
        PBar = ProgressBar('Area Perturbation Progress:', Steps = (NumberFrames-1), BarLen = 20, UpdateFreq = 1.)
        
        for i in range(NumberFrames):
            
            PBar.Update(i)
            
            xyz = traj[i].xyz[0]
            xyz_red = self.ReduceCoordinates(xyz,periodicboxvectors) # get reduced coordinates

            time_start = time.time() # time how long it takes
            
            if ReRunRefState:
                simulation0.context.setPeriodicBoxVectors(box[0],box[1],box[2])
                simulation0.context.setPositions(traj[i].xyz[0])
                simulation0.context.applyConstraints(1e-12)
                state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
                temp_PotEne_0 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            else: 
                temp_PotEne_0 = PotEne_State0_List[i]            
            
            PotEne_State0.append(temp_PotEne_0)
            
            # Find reference positions
            natoms = traj.topology.n_atoms
            reference = np.zeros([natoms,3])
            for res in traj.topology.residues:
                atomids = [atom.index for atom in res.atoms]
                com = np.mean(xyz[atomids,:], 0 )
                reference[atomids,:] = com[None,:]
            
            #calculate larger box
            if dim_max == 2:
                newpos = xyz + reference*np.array( [scaleXY1-1.0, scaleXY1-1.0, scaleZ1-1.0] )    
            elif dim_max == 1:
                newpos = xyz + reference*np.array( [scaleXY1-1.0, scaleZ1-1.0, scaleXY1-1.0] )    
            elif dim_max == 0:
                newpos = xyz + reference*np.array( [scaleZ1-1.0, scaleXY1-1.0, scaleXY1-1.0] )                   
                
            # DO EXPANSION
            #xyz_exp = self.UnReduceCoordinates(xyz_red,PBA_expansion/nanometer)
            xyz_exp = newpos
            simulation0.context.setPeriodicBoxVectors(PBA_expansion[0],PBA_expansion[1],PBA_expansion[2])
            simulation0.context.setPositions(xyz_exp)
            #simulation0.context.applyConstraints(1e-12)
            state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
            temp_PotEne_1 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            
            dU_0to1 = (temp_PotEne_1-temp_PotEne_0) 
            
            PotEne_Data_Expansion.append([temp_PotEne_0, temp_PotEne_1, dU_0to1])
            flog = open('ReRun_PotEne_State0to1_AreaExpansion.log','a')
            flog.write('{},  {},  {},  {}\n'.format(i,temp_PotEne_0, temp_PotEne_1, dU_0to1))
            flog.close()    
        
            #calculate smaller box
            if dim_max == 2:
                newpos = xyz + reference*np.array( [scaleXY2-1.0, scaleXY2-1.0, scaleZ2-1.0] )
            elif dim_max == 1:
                newpos = xyz + reference*np.array( [scaleXY2-1.0, scaleZ2-1.0, scaleXY2-1.0] )
            elif dim_max == 0:
                newpos = xyz + reference*np.array( [scaleZ2-1.0, scaleXY2-1.0, scaleXY2-1.0] )
        
            # DO CONTRACTION
            #xyz_cont = self.UnReduceCoordinates(xyz_red,PBA_contraction/nanometer)
            xyz_cont = newpos
            simulation0.context.setPeriodicBoxVectors(PBA_contraction[0],PBA_contraction[1],PBA_contraction[2])
            simulation0.context.setPositions(xyz_cont)
            #simulation0.context.applyConstraints(1e-12)
            state = simulation0.context.getState(getEnergy=True,enforcePeriodicBox=False)
            temp_PotEne_2 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            
            dU_1to0 = (temp_PotEne_2-temp_PotEne_0) 
            
            PotEne_Data_Contraction.append([temp_PotEne_0, temp_PotEne_2, dU_1to0])
            flog = open('ReRun_PotEne_State0to1_AreaContraction.log','a')
            flog.write('{},  {},  {},  {}\n'.format(i,temp_PotEne_0, temp_PotEne_2, dU_1to0))
            flog.close()    
        
        PBar.Clear()
    
        return PotEne_Data_Expansion, PotEne_Data_Contraction, PotEne_State0
    
    
    def Reweight(self):
        '''
            Reweight generates the dUs = U_1-U_0 where 1 and 0 denote the two states, 0 being the reference state.
        '''
        
        if len(self.traj_list) == 2: 
        #decide if we are doing single-state or two-state reweighting
            self.number_states = 2
                
        for traj_index, traj in enumerate(self.traj_list):
            
            ''' CHEMICAL POTENTIAL CALCULATION '''
            if self.derivative == 'particle': # calculating chemical potential 
                print('Applying a particle perturbation...')
                
                if traj_index == 0: # do insertions first 
                    # State0to1 is going from high to low entropy (i.e. inserting a molecule)
                    if self.NumberInsertionsPerFrame > 0: 
                        # Generate Molecular Conformation Library
                        NewPositions_Lists = self.GenerateMoleculeConfigsInsertions(traj,self.NumberDrawsPerFrame,self.NumberInsertionsPerFrame)                
                    
                        # Do Insertions
                        PotEne_Data_0to1_List, PotEne_State0_List = self.PerformInsertions(traj,self.states_list[0],self.states_list[1],NewPositions_Lists,self.thermo_files_list[0],self.ReRunRefState)
                    
                        if len(PotEne_State0_List) > 0 and len(self.thermo_files_list[0]) != 0:
                            self.thermo_files_list.append(PotEne_State0_List) # this way do not have to rerun ref. state twice 
                    
                        if self.number_states == 1: # Approximate state 1 --> 0 with deletions, i.e. do not have two thermo_files in the thermo_files_list
                            PotEne_Data_1to0_List = self.PerformDeletions(traj,self.states_list[0],self.states_list[2],self.NumberDeletionsPerFrame,self.thermo_files_list[0],self.ReRunRefState)
                    else: #if dont do insertion, rerun state to get PE
                        PotEne_Data_1to0_List = self.PerformDeletions(traj,self.states_list[0],self.states_list[2],self.NumberDeletionsPerFrame,self.thermo_files_list[0],True) 
                elif self.number_states == 2 and traj_index == 1: # Have actual state 1 trajectory (N+1) and need to reweight to N 
                    PotEne_Data_1to0_List = self.PerformDeletions(traj,self.states_list[1],self.states_list[0],self.NumberDeletionsPerFrame,self.thermo_files_list[1],self.ReRunRefState)
                
            ''' PRESSURE CALCULATION '''
            if self.derivative == 'volume': # calculating pressure
                print('Applying a volume perturbation...')            
                                
                if traj_index == 0: # do expansion and contraction simulatenously if only one traj
                    # State0to1 is going from high to low entropy (i.e. inserting a molecule, expanding a box)
                    
                    # Do Expansion and Contraction
                    PotEne_Data_0to1_List, PotEne_Data_1to0_List, PotEne_State0_List = self.ApplyVolumePerturbation(traj,self.states_list[0],self.thermo_files_list[0],self.ReRunRefState)
                    
                    if len(PotEne_State0_List) > 0 and len(self.thermo_files_list[0]) != 0:
                        self.thermo_files_list.append(PotEne_State0_List) # this way do not have to rerun ref. state twice 
                                    
                
                elif self.number_states == 2 and traj_index == 1: # Have actual state 1 trajectory (N+1) and need to reweight to N 
                    ''' Need to complete reweighting from two explicit states for the volume. '''
                    pass
                    
                self.dU_0to1 = PotEne_Data_0to1_List # expansion
                self.dU_1to0 = PotEne_Data_1to0_List # contraction                
                
            ''' SURFACE TENSION CALCULATION '''
            if self.derivative == 'area': # calculating the surface tension
                print('Applying an area perturbation...')    
                
                if traj_index == 0: # do expansion and contraction simulatenously if only one traj
                    # State0to1 is going from high to low entropy (i.e. inserting a molecule, expanding a box)
                    
                    # Do Expansion and Contraction
                    PotEne_Data_0to1_List, PotEne_Data_1to0_List, PotEne_State0_List = self.ApplyAreaPerturbation(traj,self.states_list[0],self.thermo_files_list[0],self.ReRunRefState)
                    
                    if len(PotEne_State0_List) > 0 and len(self.thermo_files_list[0]) != 0:
                        self.thermo_files_list.append(PotEne_State0_List) # this way do not have to rerun ref. state twice 
                                    
                
                elif self.number_states == 2 and traj_index == 1: # Have actual state 1 trajectory (N+1) and need to reweight to N 
                    ''' Need to complete reweighting from two explicit states for the volume. '''
                    pass
                    
                self.dU_0to1 = PotEne_Data_0to1_List # expansion
                self.dU_1to0 = PotEne_Data_1to0_List # contraction    
            #''' ADD NEW CALCULATIONS HERE '''        
        self.dU_0to1 = PotEne_Data_0to1_List
        self.dU_1to0 = PotEne_Data_1to0_List
