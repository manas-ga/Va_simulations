#!/bin/bash

#$ -V
#$ -cwd
#$ -N TEST
#$ -t 1
#$ -tc 1
#$ -l h_rt=24:10:00
#$ -l h_vmem=100g
#$ -pe sharedmem 2
#$ -j y
#$ -o ~/std_out/


# std_out on Eddie: ~/std_out/
# std_out on AC3: /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_Control_sims.R $Param
