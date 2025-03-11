#!/bin/bash

#$ -V
#$ -cwd
#$ -N lost_vw
#$ -t 1-300
#$ -tc 3
#$ -l h_vmem=10g
#$ -pe smp64 2
#$ -j y
#$ -o /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat lost_va_IDs.txt | awk "NR==$SGE_TASK_ID"`

Rscript lost_va_script.R $Param
