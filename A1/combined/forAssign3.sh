per

#!/bin/bash
#BSUB -q hpcintrogpu
#BSUB -n 1
#BSUB -W 15
#BSUB -J forAssign3
#BSUB -o forAssign3%J.out
#BSUB -N
#BSUB -R "rusage[mem=4GB]"
#BSUB -R "span[hosts=1]"


# load the needed compiler here

module load gcc/8.3.0

Ns="64 128 256 512 1024 2048 4096 8192"

export MFLOPS_MAX_IT=1 
export MATMULT_COMPARE=0

for N in $Ns
do
    matmult_c.gcc lib $N $N $N
done

