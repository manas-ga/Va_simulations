#!/bin/bash

#$ -V
#$ -cwd
#$ -N Simplified_sims_pool_seq
#$ -t 1-100
#$ -tc 50
#$ -l h_rt=2:30:00
#$ -l h_vmem=2G
#$ -pe sharedmem 5
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
