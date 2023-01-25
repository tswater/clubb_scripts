#!/bin/sh
#SBATCH --nodes=6
#SBATCH --tasks-per-node=1
#SBATCH --job-name="clb_rn1"
#SBATCH --output="log2.txt"
#SBATCH --exclusive
#--exclude=node5,node6,node7,node8
cd ../
#time mpiexec -n 96 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "tall_2c/" -k 2 -f 0 -w 0
#time mpiexec -n 96 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "tall_1c/" -k 1 -f 0 -w 0
time mpiexec -n 96 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "tall2_cpl/" -k 2

