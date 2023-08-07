#!/bin/bash
#SBATCH -N 25
#SBATCH --tasks-per-node 2
#SBATCH --job-name=du_reeval_caesb
#SBATCH --output=out_du/du_reeval_caesb.out
#SBATCH --error=out_du/du_reeval_caesb.err
#SBATCH --time=30:00:00
#SBATCH --mail-user=bruna.mattos@aluno.unb.br
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=8

module load openmpi3/3.1.4
module spider py3-numpy/1.15.3

START="$(date +%s)"

mpirun python du_reeval_script.py

DURATION=$[ $(date +%s) - ${START} ]

echo ${DURATION}