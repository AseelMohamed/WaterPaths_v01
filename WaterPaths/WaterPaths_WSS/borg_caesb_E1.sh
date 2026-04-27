#!/bin/bash
#SBATCH --nodes=9
#SBATCH --ntasks=216
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --threads-per-core=1
#SBATCH --mem=180G
#SBATCH --partition=defq
#SBATCH --job-name=WSS_E1_s1
#SBATCH --output=BorgTest_E1_S1.out
#SBATCH --error=BorgError_E1_S1.err
#SBATCH --time=100:00:00

module purge
module load OLD/opt/all
module load development/compilers/GCC/12.1.0
module load development/componentsDev2024_openMPI_5.0.3

export OMPI_MCA_btl_base_warn_component_unused=0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun -np 216 ./FDBsimulation -T 2 -t 2080 -r 999 -d /eejit/home/moham073/WaterPaths_WSS/  -C -1 -O rof_tables/ -e 1   -U InputFiles/utilities_rdm_reeval.csv   -W InputFiles/water_sources_rdm_reeval.csv   -P InputFiles/policies_rdm_reeval.csv  -b true -o 5000 -n 50000 -E 1
