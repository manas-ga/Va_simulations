#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_test
#$ -t 1-18
#$ -tc 18
#$ -l h_vmem=50G
#$ -e /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim/b_Interim_files/std_out/
#$ -o /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_History_sim_JARROD.R $Param

