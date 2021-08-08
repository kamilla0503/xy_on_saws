#!/bin/bash
for ((i=0 ; i < 1000; i+=100))
do
  for len in 100 300 500 750 10000
  do
  sbatch --time=1-0:0 --wrap="srun ../mcmc $len $i"
  done
done
