#!/bin/bash

#SBATCH --job-name=gen_data_altexp
#SBATCH --output=gen_data.out
#SBATCH --error=runsim.err 
#SBATCH --time=06:00:00

julia gen_data.jl -r $1 -b $2 
