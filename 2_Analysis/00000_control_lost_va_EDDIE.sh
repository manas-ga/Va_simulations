#!/bin/bash

#$ -V
#$ -cwd
#$ -N dom_new_va
#$ -t 1-30
#$ -tc 30
#$ -l h_rt=03:30:00
#$ -l h_vmem=5G
#$ -pe sharedmem 12
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript lost_va_script.R $Param
