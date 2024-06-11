#!/bin/bash
input="000_parameter_grid.txt"
x=1
while IFS= read -r line
do
  echo Simulation number $x in progress...
  x=x+1
  echo Rscript 00_History_sim_JARROD.R $line
  Rscript 00_History_sim_JARROD.R $line
done < "$input"