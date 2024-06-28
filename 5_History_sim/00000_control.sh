#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_test
#$ -t 1-8
#$ -tc 4
#$ -l mem_free=450g,s_vmem=480g,h_vmem=500g
#$ -pe smp64 25
#$ -e /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/
#$ -o /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_History_sim_JARROD.R $Param

