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
traj0 = 'traj_wrapped.lammpstrj'
at = 2 #atomtype of particle want to insert
n = 1000 #number of insertion per frame
top = 'anim0.pdb' #pdb file to open mdtraj
trajwidom = 'traj_widom'
python ~/bin/widom/getWidomTraj.py -i $traj0 -at $at -n $n -top $top -o $trajwidom
while:
do
    if [-f $trajwidom]; then
        echo "checking file $PWD/$trajwidom"
        while:
        do
            if ! ['lsof | grep $PWD/$trajwidom']
            then
                break 
    fi
    sleep 0.5
done
echo "done"
