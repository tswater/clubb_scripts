#!/bin/sh
#SBATCH --nodes=6
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_tun"
#SBATCH --output="logtune.txt"
#SBATCH --exclusive
#SBATCH --exclude=node5,node6,node7,node8
cd ../
time mpiexec -n 192 python param_iterate.py

