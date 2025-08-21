#!/bin/bash

#$ -V
#$ -cwd
#$ -N Qpoolseq_100x
#$ -t 1-100
#$ -tc 100
#$ -l h_rt=01:00:00
#$ -l h_vmem=5G
#$ -pe sharedmem 2
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
