#!/bin/bash
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -W 15
#BSUB -J bench_poisson
#BSUB -o bench_poisson_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

module load clang

Ns="10 50 100 200"
CMD=poisson_g
IT=10000
TOL=0.005
TS=15

for N in $Ns
do
	/bin/time ./$CMD $N $IT $TOL $TS 4
done
