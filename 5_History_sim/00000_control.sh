#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_sim
#$ -t 1-10
#$ -tc 3
#$ -l mem_free=250g,s_vmem=350g,h_vmem=450g
#$ -pe smp64 10
#$ -j y
#$ -o /data/obbard/Va_simulations/analyses/b_Interim_files/std_out


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_History_sim_JARROD.R $Param

