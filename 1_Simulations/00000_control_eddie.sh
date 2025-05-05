#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_sim
#$ -t 1-3
#$ -tc 3
#$ -l h_vmem=4g
#$ -pe sharedmem 2
#$ -j y
#$ -o /exports/eddie/scratch/msamant/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/b_Interim_files/std_out/
# std_out on AC3: /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_Control_sims.R $Param
