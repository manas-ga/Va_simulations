#!/bin/bash

#$ -V
#$ -cwd
#$ -N BC_V_A
#$ -t 1-150
#$ -tc 15
#$ -l mem_free=50g,s_vmem=65g,h_vmem=75g
#$ -pe smp64 9
#$ -j y
#$ -o /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat lost_va_IDs.txt | awk "NR==$SGE_TASK_ID"`

Rscript bc_va_script.R $Param
