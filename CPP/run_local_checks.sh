#!/bin/bash
for i in 0 3 6 9 12 15
do
  for len in 5 8 11
  do
  sbatch --time=1-0:0 --wrap="srun ../mcmc $len $i"
  done
done
