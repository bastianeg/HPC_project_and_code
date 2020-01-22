#!/bin/bash
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -W 15
#BSUB -J gpu4test
#BSUB -o gpu4test_%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"

# load the needed compiler here
module load cuda/10.2
module load gcc/8.3.0

ELEMPT="2,4,8,16,32"
CMD=matmult_f.nvcc
N=800
TYPE=gpu4
export MFLOPS_MAX_IT=1
export MATMULT_COMPARE=0

for ele in $ELEMPT
do
    NUM_ELEM_PER_THREAD=$ele ./$CMD $TPYE $N $N $N
done
