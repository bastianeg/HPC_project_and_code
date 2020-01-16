#!/bin/bash
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -W 15
#BSUB -J bench_poisson
#BSUB -o bench_poisson_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here (uncomment and adjust compiler and version!)


module load gcc

N="100"
CMD=poisson_j
IT=2000
TOL=0.05
TS=15


./$CMD $N $IT $TOL $TS 4

