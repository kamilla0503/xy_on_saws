#!/bin/bash
for ((i=1300 ; i < 1400; i+=5))
do
  for len in 2000 3000 5000
  do
  sbatch --time=10-0:0 --wrap="srun ../mcmc $len $i"
  done
done
for ((i=1300 ; i < 1400; i+=5))
do
  for len in 1000
  do
  sbatch --time=1-0:0 --wrap="srun ../mcmc $len $i"
  done
done
