#!/bin/bash

#$ -V
#$ -cwd
#$ -N TEST_dom
#$ -t 1-30
#$ -tc 30
#$ -l h_rt=5:00:00
#$ -l h_vmem=25G
#$ -pe sharedmem 12
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
