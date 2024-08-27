#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_sim
#$ -t 1-27
#$ -tc 9
#$ -l mem_free=250g,s_vmem=350g,h_vmem=450g,h=!bigyin
#$ -pe smp64 1
#$ -j y
#$ -o /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


# std_out on Eddie: /exports/eddie/scratch/msamant/Va_simulations/b_Interim_files/std_out/
# std_out on AC3: /data/obbard/Va_simulations/analyses/b_Interim_files/std_out/


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_History_sim_JARROD.R $Param
