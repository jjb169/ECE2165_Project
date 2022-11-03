#!/bin/bash
module purge
module load gcc/9.2.0
module load openmpi/4.0.5

mpic++ bc_mpi.cpp -o bc_mpi -O3

