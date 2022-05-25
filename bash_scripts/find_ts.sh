#!/bin/sh
#SBATCH --nodes=3
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_ru1"
#SBATCH --output="logc1.txt"
#SBATCH --exclusive
cd ../
mpiexec -n 92 python sgp_iterate.py -k 1 -l 400 -w 0 -f 0 -y 30 -d "sgp_1c_dt30/"
mpiexec -n 92 python sgp_iterate.py -k 1 -l 400 -w 0 -f 0 -y 12 -d "sgp_1c_dt12/"
mpiexec -n 92 python sgp_iterate.py -k 1 -l 400 -w 0 -f 0 -y 6 -d "sgp_1c_dt6/"
mpiexec -n 92 python sgp_iterate.py -k 1 -l 400 -w 0 -f 0 -y 1 -d "sgp_1c_dt1/"

