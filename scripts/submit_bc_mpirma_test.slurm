#!/bin/bash
#SBATCH --nodes=4
#SBATCH --cluster=mpi
#SBATCH --partition=mpi
#SBATCH --time=8:00:00
#SBATCH --job-name=bc_mpi
#SBATCH --output=outputs/mpi-CA-GrQc-%j.txt
#SBATCH --mail-type=BEGIN,END,FAIL

# Navigate to correct directory
cd ~/ece2165

# Purge and Load modules
module purge
module load gcc/9.2.0
module load openmpi/4.0.5

# Compile program
./scripts/compile_mpi.sh

# 1*1
# mpirun -n 1 --mca btl_openib_allow_ib 1 ./bc_mpi 2 CA-GrQc
# 2*1
# mpirun -n 2 --mca btl_openib_allow_ib 1 ./bc_mpi 2 CA-GrQc
# 4*1
mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 2 CA-GrQc

crc-job-stats.py
