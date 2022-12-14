#!/bin/bash
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=1
#SBATCH --cluster=mpi
#SBATCH --partition=mpi
#SBATCH --time=2:00:00
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

# 1E-6 ERROR RATE
echo "<><><> 1e-6 Error Rate <><><>"
mpirun -n 2 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.000001
mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.000001
mpirun -n 8 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.000001
mpirun -n 16 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.000001
mpirun -n 32 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.000001

# 1E-5 ERROR RATE
echo "<><><> 1e-5 Error Rate <><><>"
mpirun -n 2 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.00001
mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.00001
mpirun -n 8 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.00001
mpirun -n 16 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.00001
mpirun -n 32 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.00001

# 1E-4 ERROR RATE
echo "<><><> 1e-4 Error Rate <><><>"
mpirun -n 2 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.0001
mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.0001
mpirun -n 8 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.0001
mpirun -n 16 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.0001
mpirun -n 32 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.0001

# 1E-3 ERROR RATE
echo "<><><> 1e-3 Error Rate <><><>"
mpirun -n 2 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.001
mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.001
mpirun -n 8 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.001
mpirun -n 16 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.001
mpirun -n 32 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.001

# 1E-2 ERROR RATE
echo "<><><> 1e-2 Error Rate <><><>"
mpirun -n 2 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.01
mpirun -n 4 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.01
mpirun -n 8 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.01
mpirun -n 16 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.01
mpirun -n 32 --mca btl_openib_allow_ib 1 ./bc_mpi 5 CA-GrQc 0.01

crc-job-stats.py