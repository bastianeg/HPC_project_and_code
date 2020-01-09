#!/bin/bash
#BSUB -J make_blk_plots
#BSUB -o make_blk_plots%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 15

PERM="blk"
SMALLSIZES="16 32 64 128"
BIGSIZES="256 512 1024 2048"

EXECUTABLE=matmult_c.gcc
BLKSIZE=32
HWCOUNT="-h dch,on,dcm,on,l2h,on,l2m,on"


export MATMULT_COMPARE=0


export MFLOPS_MAX_IT=500
for size in $SMALLSIZES
do
    ./$EXECUTABLE $PERM $size $size $size $BLKSIZE
done

export MFLOPS_MAX_IT=1
for size in $BIGSIZES
do
    ./$EXECUTABLE $PERM $size $size $size $BLKSIZE
done
exit 0