#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_test
#$ -t 1-10
#$ -tc 10
#$ -l mem_free=250G
#$ -l s_vmem=350G
#$ -l h_vmem=400G
#$ -pe smp 64
#$ -e /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/
#$ -o /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_History_sim_JARROD.R $Param

