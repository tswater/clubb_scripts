#!/bin/sh
#SBATCH --nodes=2
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_rn1"
#SBATCH --output="log2.txt"
#SBATCH --exclusive
#--exclude=node5,node6,node7,node8
cd ../
#time mpiexec -n 64 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "sgp_2c_p/" -k 2 -f 0 -w 0
#time mpiexec -n 64 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "sgp_1c_p/" -k 1 -f 0 -w 0
time mpiexec -n 64 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "sgp_cpl_pu2/" -k 2

