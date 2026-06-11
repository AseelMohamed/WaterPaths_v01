#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks=252
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem=0
#SBATCH --partition=rome
#SBATCH --job-name=MinMax
#SBATCH --output=MinMax_%j.out
#SBATCH --error=MinMax_%j.err
#SBATCH --time=20:00:00

module purge
module load 2025
module load OpenMPI/5.0.8-GCC-14.3.0

unset PRTE_MCA_hwloc_default_binding_policy

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mpirun -np $SLURM_NTASKS \
  --map-by ppr:${SLURM_NTASKS_PER_NODE}:node:PE=${SLURM_CPUS_PER_TASK} \
  --bind-to core \
  ./FDBsimulation -T 1 -t 2086 -r 999 \
  -d /gpfs/scratch1/shared/amohamed/WaterPaths_WSS/ \
  -C -1 -O rof_tables/ -b true \
  -U InputFiles/utilities_rdm.csv \
  -W InputFiles/water_sources_rdm.csv \
  -P InputFiles/policies_rdm.csv \
  -e 0 -o 1000 -n 50000 -E 6 -w 1