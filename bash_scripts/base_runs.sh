#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_66"
#SBATCH --output="log_keep.txt"
#SBATCH --exclusive
#--exclude=node5,node6,node7,node8
cd ../
#time mpiexec -n 92 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "mar27fit_2c/" -k 2 -f 0 -w 0
#time mpiexec -n 92 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "mar27fit_1c/" -k 1 -f 0 -w 0
#time mpiexec -n 92 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/sgp_iterate.py -d "mar27fit_cpl/" -k 2
time mpiexec -n 32 python /home/tsw35/soteria/clubb/clubb_scripts/iterators/les_fit_iter.py -d "apr25_dz60dt06" -z 60 -t 6
