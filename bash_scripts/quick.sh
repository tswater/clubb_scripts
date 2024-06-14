#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name="clb_qik"
#SBATCH --output="log.txt"
#SBATCH --exclusive
#SBATCH --exclude=node5,node6,node7,node8

mpiexec -n 16 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/quick.py

