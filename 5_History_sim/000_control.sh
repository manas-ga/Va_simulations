#!/bin/bash

input="000_parameter_grid.txt"

while IFS= read -r line
do
  param=$line
  echo $param
  echo abc
 # Rscript 00_History_sim_JARROD.R $param
done < "$input"
