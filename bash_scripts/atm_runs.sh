#!/bin/sh
#SBATCH --nodes=3
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_atm"
#SBATCH --output="logatm.txt"
#SBATCH --exclusive
#SBATCH --exclude=node5,node6,node7,node8
cd ../
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_1c/" -k 1 -f 0 -w 0 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_nocpl/" -k 2 -f 0 -w 0 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.00001/" -t 0.00001 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.00005/" -t 0.00005 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.0001/" -t 0.0001 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.00015/" -t 0.00015 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.0002/" -t 0.0002 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.00025/" -t 0.00025 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.0003/" -t 0.0003 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.00035/" -t 0.00035 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.0004/" -t 0.0004 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.00045/" -t 0.00045 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
time mpiexec -n 92 python sgp_iterate.py -d "const_atm/sgp_cpl_0.0005/" -t 0.0005 -k 2 -a "/home/tsw35/soteria/clubb/data/arm_model_20170721.in"
