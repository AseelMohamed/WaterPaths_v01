#!/bin/bash
#SBATCH -n 1 -N 1 
#SBATCH --job-name=callgrind
#SBATCH --output=output/vallgrind.out
#SBATCH --error=output/vallgrind.err
#SBATCH --time=24:00:00
#SBATCH --exclusive
module load valgrind

export OMP_NUM_THREADS=16                       # tell the cube to use all 16 core within the node to run this
cd $SLURM_SUBMIT_DIR

# function call below, note that the C flag is 1 (change the -d flag to be your directory)
valgrind --tool=memcheck --leak-check=full --track-origins=yes --log-file=valgrind_memcheck.out ./triangleSimulation\
	-T 16\
        -t 2080\
        -r 1\
        -d /scratch/spec1058/WaterPaths/\
        -C 1\
	-m 0\
	-s sample_solutions.csv\
        -O rof_tables_valgrind/\
        -e 0\
        -U TestFiles/utilities_rdm.csv\
        -W TestFiles/water_sources_rdm.csv\
        -P TestFiles/policies_rdm.csv\
	-p false


