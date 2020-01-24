#!/bin/bash
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -W 15
#BSUB -J gpu4test
#BSUB -o gpu4test.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"


module load cuda/10.2
module load gcc/8.3.0

ELEMPT="2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32"
BS="2 4 8 16 32"
CMD=matmult_f.nvcc
N=1024

export MFLOPS_MAX_IT=1
export MATMULT_COMPARE=0

for b in $BS
do
BLOCK_SIZE=$b ./$CMD gpu2 $N $N $N
for ele in $ELEMPT
do
    NUM_ELEM_PER_THREAD=$ele BLOCK_SIZE=$b ./$CMD gpu4 $N $N $N
done
done
