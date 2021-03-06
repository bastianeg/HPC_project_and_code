#!/bin/bash
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -W 15
#BSUB -J bench_poisson
#BSUB -o bench_poisson_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here (uncomment and adjust compiler and version!)
module load clang/9.0.0

THREADS="1 2 4 8 12 16 20 24"
CMD=poisson_jp
N=100
IT=2000
TOL=0.05
TS=15

for t in $THREADS
do
    OMP_NUM_THREADS=$t /bin/time ./$CMD $N $IT $TOL $TS 
done
