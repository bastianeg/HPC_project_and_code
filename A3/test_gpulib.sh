#!/bin/bash
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -W 15
#BSUB -J gpulibtest
#BSUB -o gpulibtest_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here
module load cuda/10.2
module load gcc/8.3.0

N="64 128 256 512 1024"
CMD=matmult_f.nvcc
TYPE=gpulib
export MFLOPS_MAX_IT=1
export MATMULT_COMPARE=0

for n in $N
do
    ./$CMD $TPYE $n $n $n
done
