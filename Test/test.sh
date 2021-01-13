#!/bin/bash

#SBATCH --time=36:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=9
#SBATCH --mem-per-cpu=9000M   # memory per CPU core

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load g16


~/TSS/bin/NEWTSS ~/compute/TSS/Beta/Methane/Pt_Methane_MeCN_MeCN_PMe3/one.in ~/compute/TSS/Beta/Methane/Pt_Methane_MeCN_MeCN_PMe3/Pt_Methane_MeCN_MeCN_PMe3.xyz