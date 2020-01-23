#!/bin/bash
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 15
#BSUB -J gpu2test
#BSUB -o gpu2test_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"


# load the needed compiler here
module load cuda/10.2
module load gcc/8.3.0

Ns="64 128 256 512 1024"
CMD=matmult_f.nvcc
export MFLOPS_MAX_IT=1 
export MATMULT_COMPARE=0

for N in $Ns
do
    ./$CMD gpu2 $N $N $N
done

