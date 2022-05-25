#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_rn0"
#SBATCH --output="log0.txt"
#SBATCH --exclusive
#SBATCH --exclude=node5,node6,node7,node8
cd ../
time mpiexec -n 32 python sgp_iterate.py -d "sgp_1c/" -k 1 -f 0 -w 0
time mpiexec -n 32 python sgp_iterate.py -d "sgp_nocpl_m/" -k 2 -f 0 -w 0
time mpiexec -n 32 python sgp_iterate.py -d "sgp_cpl_m/" -k 2

