#!/bin/bash

#$ -V
#$ -cwd
#$ -N poolseq_simpl_noQ_500x
#$ -t 1-100
#$ -tc 100
#$ -l h_rt=1:30:00
#$ -l h_vmem=5G
#$ -pe sharedmem 2
#$ -j y
#$ -o ~/std_out/

Param=`cat analysis_param_grid_eddie.txt | awk "NR==$SGE_TASK_ID"`

Rscript Analysis_script.R $Param
