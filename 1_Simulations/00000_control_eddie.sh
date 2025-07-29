#!/bin/bash

#$ -V
#$ -cwd
#$ -N Set_27
#$ -t 1-100
#$ -tc 100
#$ -l h_rt=01:20:00
#$ -l h_vmem=10g
#$ -pe sharedmem 2
#$ -j y
#$ -o ~/std_out/


# std_out on Eddie: ~/std_out/
# std_out on AC3: /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_Control_sims.R $Param
