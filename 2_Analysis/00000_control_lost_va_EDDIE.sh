#!/bin/bash

#$ -V
#$ -cwd
#$ -N Set_23_new_va
#$ -t 1-300
#$ -tc 100
#$ -l h_rt=00:12:00
#$ -l h_vmem=5G
#$ -pe sharedmem 2
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript lost_va_script.R $Param
