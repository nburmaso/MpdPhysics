#!/bin/sh
# submit a desired number of jobs consisting of several subjobs
# note: tune sbatch according to your cluster CPUs for optimal usage,
#       e.g. `sbatch -n 28 srun.sh` if a CPU has 28 threads
for (( I=1; I<100; I++ ))
  do export I
  echo starting batch job $I
  sbatch -n 28 srun.sh
done
