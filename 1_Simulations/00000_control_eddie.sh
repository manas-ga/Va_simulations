#!/bin/bash

#$ -V
#$ -cwd
#$ -N Set_22
#$ -t 1-100
#$ -tc 50
#$ -l h_rt=2:00:00
#$ -l h_vmem=24g
#$ -pe sharedmem 12
#$ -j y
#$ -o ~/std_out/


# std_out on Eddie: ~/std_out/
# std_out on AC3: /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_Control_sims.R $Param
