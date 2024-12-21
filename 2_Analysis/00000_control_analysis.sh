#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_al_Set_8
#$ -t 1-50
#$ -tc 8
#$ -l mem_free=250g,s_vmem=350g,h_vmem=450g,h=!bigyin
#$ -pe smp64 15
#$ -j y
#$ -o /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat analysis_param_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
