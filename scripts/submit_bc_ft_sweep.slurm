#!/bin/bash
#SBATCH --nodes=4
#SBATCH --cluster=mpi
#SBATCH --partition=mpi
#SBATCH --time=0:10:00
#SBATCH --job-name=ftmpi-bc
#SBATCH --output=outputs/ftmpi-CA-GrQc-%j.txt
#SBATCH --account=ageorge

# Navigate to correct directory
cd ~/ece2165

# Purge and Load modules
module purge
module load gcc/9.2.0
module load openmpi/4.0.5

# Compile program
mpic++ bc_mpi.cpp -o bc_mpi -O3

# 1*1
# mpirun -n 1 --mca btl_openib_allow_ib 1 ./bc_mpi 2 CA-GrQc
# 2*1
# mpirun -n 2 --mca btl_openib_allow_ib 1 ./bc_mpi 2 CA-GrQc
# 4*1
mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 2 CA-GrQc 0.001
#mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 2 small_example 0.001

crc-job-stats.py