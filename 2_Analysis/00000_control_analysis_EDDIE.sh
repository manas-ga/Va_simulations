#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_al_Set_10
#$ -t 1-100
#$ -tc 10
#$ -l h_rt=10:05:00
#$ -l h_vmem=750G
#$ -pe sharedmem 50
#$ -j y
#$ -o /exports/eddie/scratch/msamant/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat analysis_param_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
