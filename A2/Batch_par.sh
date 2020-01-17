#!/bin/bash
#BSUB -q hpcintro
#BSUB -n 8
#BSUB -W 15
#BSUB -J bench_poisson
#BSUB -o bench_poisson_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here (uncomment and adjust compiler and version!)
module load clang/9.0.0

THREADS="8"
Ns="10 50 100 200"
CMD=poisson_gsp
IT=1000
TOL=1.5e-3
TS=15

for N in $Ns
do
    OMP_NUM_THREADS=$THREADS OMP_PROC_BIND=close OMP_PLACES=cores /bin/time ./$CMD $N $IT $TOL $TS 
done
