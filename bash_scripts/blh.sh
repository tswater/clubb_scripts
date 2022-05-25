#!/bin/sh
#SBATCH --nodes=3
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_blh"
#SBATCH --output="logblh.txt"
#SBATCH --exclusive
#SBATCH --exclude=node5,node6,node7,node8
cd ../
time mpiexec -n 96 python add_blh.py

