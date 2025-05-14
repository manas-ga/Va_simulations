#!/bin/bash

#$ -V
#$ -cwd
#$ -N Set_5_re
#$ -t 1-900
#$ -tc 15
#$ -l h_rt=0:5:00
#$ -l h_vmem=2G
#$ -pe sharedmem 12
#$ -j y
#$ -o /exports/eddie/scratch/msamant/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
