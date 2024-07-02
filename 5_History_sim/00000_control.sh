#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_sim
#$ -t 1-90
#$ -tc 10
#$ -l mem_free=250g,s_vmem=350g,h_vmem=450g
#$ -pe smp64 22
#$ -e /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/
#$ -o /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_History_sim_JARROD.R $Param

