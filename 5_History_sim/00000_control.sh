#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_sim
#$ -t 1-900
#$ -tc 10
#$ -l mem_free=130g,s_vmem=120g,h_vmem=150g
#$ -pe smp64 20
#$ -e /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/
#$ -o /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Param=`cat 000_parameter_grid.txt | awk "NR==$SGE_TASK_ID"`

Rscript 00_History_sim_JARROD.R $Param

