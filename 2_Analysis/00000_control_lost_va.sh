#!/bin/bash

#$ -V
#$ -cwd
#$ -N New_V_A_set_N1
#$ -t 157-168
#$ -tc 15
#$ -l mem_free=50g,s_vmem=65g,h_vmem=75g
#$ -pe smp64 9
#$ -j y
#$ -o /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat lost_va_IDs.txt | awk "NR==$SGE_TASK_ID"`

Rscript lost_va_script.R $Param
