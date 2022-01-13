#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_ru1"
#SBATCH --output="logc1.txt"
#SBATCH --exclusive

time python cpl_ant.py
