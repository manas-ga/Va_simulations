#!/bin/bash

#$ -V
#$ -cwd
#$ -N Set_N1_re
#$ -t 1-400
#$ -tc 40
#$ -l h_rt=4:00:00
#$ -l h_vmem=22G
#$ -pe sharedmem 12
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
