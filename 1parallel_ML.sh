#!/bin/bash

#PBS -N te_pa
#PBS -l select=1:ncpus=4
#PBS -l walltime=20:00:00
#PBS -o /dev/null
#PBS -e /dev/null

#Loads R
module load R

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
mpirun -np 1 --hostfile $PBS_NODEFILE -x PATH -x LD_LIBRARY_PATH  R --no-save < 1parallel_ML.R > output 2>&1
