#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name="clb_tun"
#SBATCH --output="logtune.txt"
#SBATCH --exclusive
# --exclude=node5,node6,node7,node8
cd ../
time mpiexec -n 32 python /home/tsw35/soteria/clubb/clubb_scripts/jupyter_notebooks/les_tune_data.py
#python /home/tsw35/soteria/clubb/clubb_scripts/add_blh.py
