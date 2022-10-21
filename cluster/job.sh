#!/bin/bash

#SBATCH --partition=full

#SBATCH --job-name=IMT2112

#SBATCH --output=log.out

#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1

mpic++ main.cpp -std=c++11
time mpirun ./a.out 2000