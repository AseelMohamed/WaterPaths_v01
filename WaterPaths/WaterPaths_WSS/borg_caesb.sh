#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=2
#SBATCH --threads-per-core=1
#SBATCH --mem=0
#SBATCH --partition=rome
#SBATCH --job-name=BorgWSS2
#SBATCH --output=BorgTest5.out
#SBATCH --error=BorgError5.err
#SBATCH --time=05:00:00

module purge
module load 2025
module load OpenMPI/5.0.8-GCC-14.3.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores 

mpirun -np 36 ./FDBsimulation -T 2 -t 2086 -r 999 -d /gpfs/scratch1/shared/amohamed/WaterPaths_WSS/ \
  -C -1 -O rof_tables/ -b true \
  -U InputFiles/utilities_rdm_reeval.csv \
  -W InputFiles/water_sources_rdm_reeval.csv \
  -P InputFiles/policies_rdm_reeval.csv \
  -e 5 -o 25 -n 100 -E 2