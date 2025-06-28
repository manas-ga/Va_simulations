#!/bin/bash

#$ -V
#$ -cwd
#$ -N Full_sims_pool_seq
#$ -t 88
#$ -tc 1
#$ -l h_rt=8:00:00
#$ -l h_vmem=24G
#$ -pe sharedmem 12
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
