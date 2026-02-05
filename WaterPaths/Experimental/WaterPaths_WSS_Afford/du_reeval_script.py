# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 10:49:25 2021

@author: lbl59
"""

from mpi4py import MPI
import numpy as np
import subprocess, sys, time
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

N_RDMs = 1000 #Bruna: I think as a test 2 would be enough # DAVE- DEPENDS ON HOW MANY NODES, WE NEED MORE TOTAL RDM THAN TOTAL TASKS

OMP_NUM_THREADS = 8  #Bruna: To test only put 2 # DAVE- LETS USE 8 SINCE THAT's HOW MANY WE'LL USE FOR THE FULL RUN
N_REALIZATIONS = 1000 #Bruna: put only 10 realizations as a test, but the right amount would be realizations=1000, right?
DATA_DIR = "/scratch/spec1058/WaterPaths_jan22/" #Bruna: changed to overall WaterPaths_jan22 folder

SOLS_FILE_NAME = "./DVS_refset.csv" #Bruna: didn't put any specific DU directory,is this right?

N_NODES = 25 #Bruna: put only 2 nodes
N_TASKS_PER_NODE = 2 # Bruna: put 5 tasks per node to give 10 realizations, check if this logic is correct # DAVE - THIS SHOULD BE THE TOTAL NUMBER OF CORE ON A NODE (16) DIVIDED BY OMP_NUM_THREADS (8)
N_TASKS = int(N_TASKS_PER_NODE * N_NODES) # should be 200 # DAVE - SHOULD BE 4 
N_RDMS_PER_TASK = int(N_RDMs/N_TASKS)  # should be 5 # DAVE - SHOULD BE 2
 
SOL_NUM = 91   #Bruna: caesb policies

for i in range(N_RDMS_PER_TASK):
    current_RDM = rank + (N_TASKS * i)

#changed a lot on the commands below
    command_run_rdm = "./FDBsimulation /scratch/spec1058/WaterPaths_jan22 -T {} -t 2080 -r {} -d {} -C -1 -O rof_tables_no112/ -e 0 \
        -U InputFiles_no112/utilities_rdm_reeval.csv \
        -W InputFiles_no112/water_sources_rdm_reeval.csv \
        -P InputFiles_no112/policies_rdm_reeval.csv \
        -R {} -s {} -f 0 -l {} -p false -c false".format(OMP_NUM_THREADS, N_REALIZATIONS, DATA_DIR, current_RDM, SOLS_FILE_NAME, SOL_NUM)
    
    print(command_run_rdm)
    os.system(command_run_rdm)

comm.Barrier()