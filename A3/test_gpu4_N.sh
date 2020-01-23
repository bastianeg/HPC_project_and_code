#!/bin/bash
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 15
#BSUB -J gpu4test_N
#BSUB -o gpu4test_N_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here
module load cuda/10.2
module load gcc/8.3.0

Ns="64 128 256 512 1024"
CMD=matmult_f.nvcc
TYPE=gpu4
export MFLOPS_MAX_IT=1
export MATMULT_COMPARE=0

for N in $Ns
do
    NUM_ELEM_PER_THREAD=2 ./$CMD gpu4 $N $N $N
done
