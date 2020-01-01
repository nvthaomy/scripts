#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=100:00:00
#PBS -V
#PBS -N widom
#PBS -M my@ucsb.edu
#PBS -m ae

cd $PBS_O_WORKDIR
export PATH=/home/mnguyen/bin/lammps/lammps-12Dec18/bin/:$PATH
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_git:$PYTHONPATH
Ntot=6000
NP=4
NS=5856
N=36
L=40
inputF='polymer.in'
dcd='traj.dcd'
traj0='traj_wrapped.lammpstrj'
at=3 #atomtype of particle want to insert
n=1000 #number of insertion per frame
top='anim0.pdb' #pdb file to open mdtraj
dat0='polymer0.data'
trajwidom='traj_widom'
traj1=$trajwidom'.lammpstrj'
dat1name='Ntot'$((Ntot+1))
dat1=$dat1name'.data'
python ~/bin/widom/getWidomTraj.py -i $traj0 -at $at -n $n -top $top -o $trajwidom -dcd $dcd
python ~/bin/BSpolymerExpSolventInit.py -np $NP -ns $((NS+1)) -N $N -L $L -dat $dat1name
python ~/bin/widom/writeInput.py -i $inputF -dat $dat0 -trj $traj0 -o rerun0 -txt u0
python ~/bin/widom/writeInput.py -i $inputF -dat $dat1 -trj $traj1 -o rerun1 -txt u1
lmp_omp -sf omp -pk omp 2 -in rerun0.in -log rerun0.log
lmp_omp -sf omp -pk omp 2 -in rerun1.in -log rerun1.log
python ~/bin/widom/getChemPot.py -u0 u0.txt -u1 u1.txt -n $n


