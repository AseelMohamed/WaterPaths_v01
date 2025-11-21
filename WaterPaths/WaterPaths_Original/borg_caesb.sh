#!/bin/bash
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --nodes=2
#SBATCH --partition=defq
#SBATCH --job-name=BorgTest2_s0
#SBATCH --output=test50k_s0.out
#SBATCH --error=error50k_s0.err
#SBATCH --time=15:00:00

module purge
module load OLD/opt/all
module load development/compilers/GCC/12.1.0
module load development/componentsDev2024_openMPI_5.0.3

mpirun -np 24 ./FDBsimulation\
        -T 2\
        -t 2080\
        -r 100\
        -d /eejit/home/moham073/WaterPaths_Original/\
        -C -1\
        -O rof_tables/\
        -e 0\
        -U InputFiles/utilities_rdm_reeval.csv\
        -W InputFiles/water_sources_rdm_reeval.csv\
        -P InputFiles/policies_rdm_reeval.csv\
        -b true\
        -o 5000\
        -n 50000