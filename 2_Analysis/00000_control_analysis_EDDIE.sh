#!/bin/bash

#$ -V
#$ -cwd
#$ -N Set_14_re
#$ -t 1-100
#$ -tc 20
#$ -l h_rt=7:00:00
#$ -l h_vmem=22G
#$ -pe sharedmem 12
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
