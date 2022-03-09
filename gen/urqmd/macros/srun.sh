#!/bin/sh
#SBATCH -D path/to/your/logs
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -t 4:00:00
#SBATCH -p compute
#SBATCH -J urqmd_run
date
srun urqmd.sh
date
