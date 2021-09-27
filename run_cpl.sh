#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_runs"
#SBATCH --output="logclb.txt"
#SBATCH --exclusive

time mpiexec -n 12 python lhet_iterate.py
#time mpiexec -n 22 python core_sgp_cpl.py
