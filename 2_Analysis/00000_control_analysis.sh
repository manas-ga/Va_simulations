#!/bin/bash

#$ -V
#$ -cwd
#$ -N Set_9_rep
#$ -t 1-53
#$ -tc 8
#$ -l mem_free=250g,s_vmem=350g,h_vmem=450g,h=!bigyin
#$ -pe smp64 22
#$ -j y
#$ -o /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat analysis_param_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
