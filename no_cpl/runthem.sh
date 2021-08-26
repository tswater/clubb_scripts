#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --job-name="clb_runs"
#SBATCH --output="logclb.txt"
#SBATCH --exclusive

time mpiexec -n 32 python sgprun_0.py
time mpiexec -n 32 python sgprun_1.py
time mpiexec -n 32 python sgprun_2.py
time mpiexec -n 32 python sgprun_3.py
time mpiexec -n 32 python sgprun_4.py
time mpiexec -n 32 python sgprun_5.py
time mpiexec -n 32 python sgprun_6.py
time mpiexec -n 32 python sgprun_7.py
time mpiexec -n 32 python sgprun_8.py
time mpiexec -n 32 python sgprun_9.py
time mpiexec -n 32 python sgprun_10.py
time mpiexec -n 32 python sgprun_11.py
time mpiexec -n 32 python sgprun_12.py
time mpiexec -n 32 python sgprun_13.py
time mpiexec -n 32 python sgprun_14.py
time mpiexec -n 32 python sgprun_15.py
time mpiexec -n 32 python sgprun_16.py
time mpiexec -n 32 python sgprun_17.py
time mpiexec -n 32 python sgprun_18.py
time mpiexec -n 32 python sgprun_19.py
time mpiexec -n 32 python sgprun_20.py
time mpiexec -n 32 python sgprun_21.py
time mpiexec -n 32 python sgprun_22.py
time mpiexec -n 32 python sgprun_23.py
time mpiexec -n 32 python sgprun_24.py
time mpiexec -n 32 python sgprun_25.py
time mpiexec -n 32 python sgprun_26.py
time mpiexec -n 32 python sgprun_27.py
time mpiexec -n 32 python sgprun_28.py
time mpiexec -n 32 python sgprun_29.py
time mpiexec -n 32 python sgprun_30.py
time mpiexec -n 32 python sgprun_31.py
time mpiexec -n 32 python sgprun_32.py
time mpiexec -n 32 python sgprun_33.py
time mpiexec -n 32 python sgprun_34.py
time mpiexec -n 32 python sgprun_35.py
time mpiexec -n 32 python sgprun_36.py
time mpiexec -n 32 python sgprun_37.py
time mpiexec -n 32 python sgprun_38.py
time mpiexec -n 32 python sgprun_39.py
time mpiexec -n 32 python sgprun_40.py
time mpiexec -n 32 python sgprun_41.py
time mpiexec -n 32 python sgprun_42.py
