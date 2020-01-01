#!/bin/bash -l
#SBATCH --nodes=1 --ntasks-per-node=1
cd $SLURM_SUBMIT_DIR
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
python getRgRee.py
