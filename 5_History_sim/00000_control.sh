#!/bin/bash

#$ -V
#$ -cwd
#$ -N Vw_sim
#$ -t 1-1
#$ -tc 1
#$ -l mem_free=250g,s_vmem=350g,h_vmem=450g
#$ -pe smp64 10
#$ -j y
#$ -o /ceph/users/marun/Va_simulations/5_History_sim/b_Interim_files/std_out/


Rscript 00_History_sim_JARROD.R 1.2e-06 1.4 1.4 1000 10 3
