#!/bin/bash
#BSUB -J make_blk_plots
#BSUB -o make_blk_plots%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 15

PERM="blk"
BLKSIZES="4 8 16 24 32 36 64 128 256 320 512"
SIZE="1024"
EXECUTABLE=matmult_c.gcc

export MATMULT_COMPARE=0

export MFLOPS_MAX_IT=1
for BLKSIZE in $BLKSIZES
do
    ./$EXECUTABLE $PERM $SIZE $SIZE $SIZE $BLKSIZE
done
exit 0