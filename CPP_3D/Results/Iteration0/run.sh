#!/bin/bash
for ((i=0 ; i < 180; i+=15))
do
  for len in 50 100 200 300 
  do
  sbatch --time=1-0:0 --wrap="srun ~/saw_models/mcmc $len $i"
  done
done
