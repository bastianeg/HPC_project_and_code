#!/bin/bash
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -W 15
#BSUB -J gs
#BSUB -o gs_%J.out
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here (uncomment and adjust compiler and version!)
module load clang/9.0.0

THREADS="1 2 4 8 12 16 20 24"
N="10 20 50 100 200"
CMD=poisson_gsp
IT=2000
TOL=0.05
TS=15

for t in $THREADS
do
    for n in $N
    do
    	OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores ./$CMD $N $IT $TOL $TS 
    done
done
