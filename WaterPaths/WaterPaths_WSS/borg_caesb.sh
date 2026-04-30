#!/bin/bash
#SBATCH --nodes=10
#SBATCH --ntasks=320
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem=0
#SBATCH --partition=genoa
#SBATCH --job-name=BorgWSS
#SBATCH --output=BorgTest.out
#SBATCH --error=BorgError.err
#SBATCH --time=30:00:00

module purge
module load 2025
module load OpenMPI/5.0.8-GCC-14.3.0

export OMP_NUM_THREADS=1

mpirun -np 320 ./FDBsimulation -T 1 -t 2086 -r 999 -d /gpfs/scratch1/shared/amohamed/WaterPaths_WSS/ \
  -C -1 -O rof_tables/ -e 1 \
  -U InputFiles/utilities_rdm_reeval.csv \
  -W InputFiles/water_sources_rdm_reeval.csv \
  -P InputFiles/policies_rdm_reeval.csv \
  -b true -o 1000 -n 100000 -E 2