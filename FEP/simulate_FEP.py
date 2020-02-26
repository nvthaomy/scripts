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
import FEP_Module_HOH_gro as FEP

fname = 'nacl'
f = open("sim_{}.log".format(fname), "w")

''' State Variables ''' 
Temp = (298.15)
friction = 5*picoseconds 
f.write('temperature = {} K \n'.format(Temp))

'''  Simulation Specifications '''
useGPU = False
useLJPME = True
tail =False
nonbondedCutoff = 0.9*nanometers
nThreads = 4
#FEP vars
unit = 'kJ/mol' 
resNames = ['HOH'] #name of residues to insert in each frame
unit = 'kJ/Mol'
nDraw =  5
nInsert = 5
nDelete = 5
SetReRunRefState = False

CalcChemPot = True
CalcPressure = False
CalcSurfaceTension = False
ThermoSlice =  10
TrajSlice = 2
Warmup = 100
PotEne_Data = np.loadtxt('logNPT.txt',delimiter='\t')[::ThermoSlice,1]
traj_file = 'trajNPT.dcd' # N-particle state
top_file = '6616opc_122nacl.parm7'

f.write('use GPU {}. \n'.format(useGPU))
f.write('use LJPME {}. \n'.format(useLJPME))
f.write('use dispersion correction {}. \n'.format(tail))

''' Platform specifications. '''
# Change the argument below to "CPU" or "CUDA" or "CUDA"to specify a platform 
if useGPU:
    platform0 = Platform.getPlatformByName('CUDA') 
    properties0 = {'DeviceIndex': '0', 'Precision': 'double'}
    platform1 = Platform.getPlatformByName('CUDA') 
    properties1 = {'DeviceIndex': '0', 'Precision': 'double'}
    platform2 = Platform.getPlatformByName('CUDA') 
    properties2 = {'DeviceIndex': '0', 'Precision': 'double'}
else:
    platform0 = Platform.getPlatformByName('CPU')
    properties0 = {'Threads': str(nThreads)}
    platform1 = Platform.getPlatformByName('CPU')
    properties1 = {'Threads': str(nThreads)}
    platform2 = Platform.getPlatformByName('CPU') 
    properties2 = {'Threads': str(nThreads)}

''' Input files and system object creation. '''

#gro = GromacsGroFile('input.gro')
top0 = AmberPrmtopFile('6616opc_122nacl.parm7') 
top1 = AmberPrmtopFile('6617opc_122nacl.parm7')
top2 = AmberPrmtopFile('6615opc_122nacl.parm7')
inpcrd = AmberInpcrdFile('6616opc_122nacl.crd')
#gro = gromacs.GromacsGroFile.parse('box.gro')
top0.box = inpcrd.boxVectors 
top1.box = inpcrd.boxVectors
top2.box = inpcrd.boxVectors
system0 = top0.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, ewaldErrorTolerance=0.001, rigidWater=True) #constraints=HBonds, 
system1 = top1.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, ewaldErrorTolerance=0.001, rigidWater=True) #constraints=HBonds, 
system2 = top2.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, ewaldErrorTolerance=0.001, rigidWater=True) #constraints=HBonds, 

if useLJPME:
    nbm = NonbondedForce.LJPME
else:
    nbm = NonbondedForce.PME

''' Get the different forces. '''
farray = [f for ii, f in enumerate(system0.getForces()) if isinstance(f,NonbondedForce)]
fnb = farray[0]
fnb.setNonbondedMethod(nbm) # set non-bonded method
f.write("Nonbonded method: {} \n".format(fnb.getNonbondedMethod()))

# Disable the long-range vdW correction
for i, frc in enumerate(system0.getForces()):
    if (isinstance(frc, NonbondedForce)):
        f.write('The dispersion correction setting is: {}.\n'.format(frc.getUseDispersionCorrection()))
        frc.setUseDispersionCorrection(True)
        f.write('The dispersion correction setting is: {}.\n'.format(frc.getUseDispersionCorrection()))


''' Start of FEP Stuff '''
traj_load = md.load(traj_file,top=top_file, stride=TrajSlice)
traj_load = traj_load[Warmup::]
PotEne_Data = PotEne_Data[Warmup::]
nFrames = len(traj_load)
print('Number of frames in trajectory: {}'.format(len(traj_load)))
print('Number of entries in thermo. file: {}'.format(len(PotEne_Data)))

if len(traj_load) != len(PotEne_Data):
    print("WARNING: The number of entries in the trajectory and thermo. file do not match!")

# State 0: N ; re-defined here just for clarity
system0 = top0.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, ewaldErrorTolerance=1E-7, rigidWater=True)
integrator0 = LangevinIntegrator(Temp*kelvin, friction, 0.002*picoseconds)
simulation0 = Simulation(top0.topology, system0, integrator0, platform0, properties0)
    
if CalcChemPot:
    # State 1: N+1
    system1 = top1.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, ewaldErrorTolerance=1E-7, rigidWater=True)
    integrator1 = LangevinIntegrator(Temp*kelvin, friction, 0.002*picoseconds)
    simulation1 = Simulation(top1.topology, system1, integrator1, platform1, properties1)

    # State 1: N-1
    system2 = top2.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, ewaldErrorTolerance=1E-7, rigidWater=True)
    integrator2 = LangevinIntegrator(Temp*kelvin, friction, 0.002*picoseconds)
    simulation2 = Simulation(top2.topology, system2, integrator2, platform2, properties2)

# Setup the FEP_Object
method = 'OB' # 'OB' == Optimal-Bennetts; currently specified, but you actually call different methods below on the FEP_Object
derivative = 'particle' # type of perturbation to perform
traj_list = [traj_load] # trajectory objects from MDTraj
thermo_files_list = [PotEne_Data] # IF you already have the potential energy data you can load in here. 

if CalcChemPot:
    # OpenMM simulation objects (might be memory intensive)
    states_list = [simulation0,simulation1,simulation2] 
if (CalcPressure and not CalcChemPot) or (CalcSurfaceTension and not CalcChemPot):
    # currently only need one simulation object for pressure estimate
    states_list = [simulation0] 

# Create an FEP object
FEP_Object = FEP.FEP(method,derivative,traj_list,states_list,thermo_files_list)

if CalcChemPot: 
    ''' Calculate the Chemical Potential (Excess) '''
    time1 = time.time()
    # Set a variaty of object attributes
    FEP_Object.SetMethod('OB') # optimal-Bennett's method
    FEP_Object.SetDerivative('particle') # particle derivative
    FEP_Object.SetResidueName(resNames) # the residue type to be inserted/deleted;
    # if SetResidueName(''), then defaults to the first residue in residue list.
    FEP_Object.SetNumberOfDrawsPerFrame(nDraw) # How many molecule samples to draw from on each frame to build insertion library.
    FEP_Object.SetNumberInsertionsPerFrame(nInsert)
    FEP_Object.SetNumberDeletionsPerFrame(nDelete)
    FEP_Object.SetReRunRefState(SetReRunRefState)
    FEP_Object.SetTemperature(Temp) # units Kelvin
    
    
    # Names and Numbers of each residue, and the number of atoms in the residue for each trajectory
    print('Names and Number of each residue in each trajectory:')
    print(FEP_Object.GetResidueNames())

    # Perform the reweighting
    FEP_Object.Reweight()
    
    print("Number of Residues used during insertions: {}".format(len(FEP_Object.ResidueIndiceList_Insertions)))
    print("Number of Residues used during deletions: {}".format(len(FEP_Object.ResidueIndiceList_Deletions)))
    
    # Calculate the free energy difference
    FEP_Object.SetBennettsConstant(-33.8) # Initial guess, KJ/mole
    FEP_Object.Optimal_Bennetts() # Run optimal-Bennetts
    time2 = time.time()
    # Report the free energy change in KJ/mole
    s = 'Excess chemical potential calculation for {}'.format(resNames)
    s += '\nNumber of insertions: {}\nNumber of deletions: {}\nNumber of traj frames: {}'.format(nInsert,nDelete,nFrames)
    s += "\ndF: %4.4f +/- %2.5f %s"%(FEP_Object.dF, FEP_Object.dF_stderr,unit)
    s += "\nBennett's Constant: %4.4f %s"%(float(FEP_Object.dF),unit)
    s += "\nTake %5.3f hours to finish the calculation."%((time2-time1)/3600)
    print(s)
    f = open('chemicalPot.dat','w')
    f.write(s)

if CalcPressure: 
    ''' Calculating the Pressure '''
    # Reset object for P-calculation
    FEP_Object.SetMethod('OB') # still using optimal-Bennett's method
    FEP_Object.SetDerivative('volume') # now using volume derivative
    FEP_Object.SetPerturbation(0.000005) # set the perturbation to apply to the volume
    FEP_Object.SetIdealGasDOFs(900) # there are 3037 water molecules in this system
    FEP_Object.SetReRunRefState(True)
    FEP_Object.SetTemperature((273+200)) # units Kelvin

    # Call varying methods on the FEP object
    FEP_Object.Reweight()

    print('Volume_State0: {0:3.8f} nm**3'.format(float(FEP_Object.Perturbation_State0)))
    print('Volume_State1: {0:3.8f} nm**3'.format(float(FEP_Object.Perturbation_State1)))
    print('Volume_State2: {0:3.8f} nm**3'.format(float(FEP_Object.Perturbation_State2)))

    FEP_Object.SetBennettsConstant(1.0) # Initial guess, KJ/mole
    FEP_Object.Optimal_Bennetts() # Run optimal-Bennetts

    # Report the free energy change in KJ/mole
    print("dF: {0:4.4f} +/- {1:2.5f} KJ/mole".format(FEP_Object.dF,FEP_Object.dF_stderr))
    print("Bennett's Constant: {0:4.4f} KJ/mole".format(float(FEP_Object.BennettsConstant)))

    FEP_Object.CalcPex() # Calculate excess pressure, run after running a FEP method (i.e. 'OB','NOB',etc) on the FEP_Object
    FEP_Object.CalcPid() # Calculate ideal pressure

    # Report the excess pressure and ideal pressure and total pressure in atm
    print("P_id: {0:4.5f} atm".format(FEP_Object.Pid))
    print("P_excess: {0:4.5f} +/- {1:4.6f} atm".format(FEP_Object.Pex,FEP_Object.Pex_stderr))
    print("P_total:  {0:4.5f} +/- {1:4.6f} atm".format((FEP_Object.Pex+FEP_Object.Pid),FEP_Object.Pex_stderr))

if CalcSurfaceTension: 
    ''' Calculating the Pressure '''
    # Reset object for P-calculation
    FEP_Object.SetMethod('OB') # still using optimal-Bennett's method
    FEP_Object.SetDerivative('area') # now using area derivative
    FEP_Object.SetPerturbation(0.0001) # set the perturbation to apply to the area
    FEP_Object.SetReRunRefState(False)
    FEP_Object.SetTemperature(298) # units Kelvin

    # Call varying methods on the FEP object
    FEP_Object.Reweight()

    print('Area_State0: {0:3.8f} nm**2'.format(float(FEP_Object.Perturbation_State0)))
    print('Area_State1: {0:3.8f} nm**2'.format(float(FEP_Object.Perturbation_State1)))
    print('Area_State2: {0:3.8f} nm**2'.format(float(FEP_Object.Perturbation_State2)))

    FEP_Object.SetBennettsConstant(1.) # Initial guess, KJ/mole
    FEP_Object.Optimal_Bennetts() # Run optimal-Bennetts

    # Report the free energy change in KJ/mole
    print("dF: {0:4.4f} +/- {1:2.5f} KJ/mole".format(FEP_Object.dF,FEP_Object.dF_stderr))
    print("Bennett's Constant: {0:4.4f} KJ/mole".format(float(FEP_Object.BennettsConstant)))

    FEP_Object.CalcSurfaceTension() # Calculate surface tension, run after running a FEP method (i.e. 'OB','NOB',etc) on the FEP_Object

    # Report the excess pressure and ideal pressure and total pressure in atm
    print("Sigma:  {0:4.5f} +/- {1:4.6f} J/m^2".format(FEP_Object.SurfaceTension,FEP_Object.SurfaceTension_stderr))
