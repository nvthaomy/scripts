#!/bin/bash
# ask for 16 cores on two nodes
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH --time=700:00:00
#SBATCH --job-name=5MNacl_FEP
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=my@ucsb.edu

module load intel/18

cd $SLURM_SUBMIT_DIR
/bin/hostname
OPENMM_CPU_THREADS=6
export PATH="/home/mnguyen/miniconda3/bin/python:$PATH"
python simulate_FEP.py 

