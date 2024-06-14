#!/bin/sh
#SBATCH --nodes=2
#SBATCH --tasks-per-node=2
#SBATCH --job-name="sfc_rep"
#SBATCH --output="log.txt"
#SBATCH --exclusive
#SBATCH --exclude=node5,node6,node7,node8
cd ../
time mpiexec -n 64 python make_sfc_5k.py
