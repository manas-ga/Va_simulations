#!/bin/bash

#$ -V
#$ -cwd
#$ -N S26_500x
#$ -t 1-100
#$ -tc 100
#$ -l h_rt=10:00:00
#$ -l h_vmem=25G
#$ -pe sharedmem 12
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
