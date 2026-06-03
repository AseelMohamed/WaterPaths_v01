#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=72
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem=0
#SBATCH --partition=rome
#SBATCH --job-name=BorgWSS
#SBATCH --output=BorgTest_%j.out
#SBATCH --error=BorgError_%j.err
#SBATCH --time=03:00:00

module purge
module load 2025
module load OpenMPI/5.0.8-GCC-14.3.0

# Not needed for pure MPI, but harmless to keep
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
  -U InputFiles/utilities_rdm_reeval.csv \
  -W InputFiles/water_sources_rdm_reeval.csv \
  -P InputFiles/policies_rdm_reeval.csv \
  -e 0 -o 500 -n 1000 -E 1