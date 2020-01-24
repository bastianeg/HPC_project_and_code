#!/bin/bash
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -W 15
#BSUB -J multipas
#BSUB -o multipas.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=2:mode=exclusive_process"


module load cuda/10.2
module load gcc/8.3.0

CMD=jacobimulti
N=10
IT=1000
TOL=1
TSTART=15

./$CMD $N $IT $TOL $TSTART

