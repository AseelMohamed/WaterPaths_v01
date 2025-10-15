#!/bin/bash
#SBATCH -n 45 -N 15
#SBATCH --job-name=caesb_test
#SBATCH --output=output/last_no112.out
#SBATCH --error=output/last_no112.err
#SBATCH --time=60:00:00
#SBATCH --exclusive

export OMP_NUM_THREADS=5
cd $SLURM_SUBMIT_DIR

time mpirun -np 45 ./FDBsimulation\
	-T 5\
        -t 2080\
        -r 999\
        -d /scratch/spec1058/WaterPaths_jan22/\
        -C -1\
        -O rof_tables_no112/\
        -e 4\
        -U InputFiles_no112/utilities_rdm.csv\
        -W InputFiles_no112/water_sources_rdm.csv\
        -P InputFiles_no112/policies_rdm.csv\
	-b true\
	-n 30000\
	-o 5000
