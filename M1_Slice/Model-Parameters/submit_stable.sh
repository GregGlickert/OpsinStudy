#!/bin/bash
#SBATCH --job-name=opto
#SBATCH --time=06:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --mem-per-cpu=2G
#SBATCH --output=../Run-Storage/run_stable.o
#SBATCH --error=../Run-Storage/run_stable.e

# Additional user-defined commands
module load mpich-x86_64-nopy

export OUTPUT_DIR=../Run-Storage/opto_stable
mpirun nrniv -mpi -python run_opto_network.py simulation_config_stable.json False
