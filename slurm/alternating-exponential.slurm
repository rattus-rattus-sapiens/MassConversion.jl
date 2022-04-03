#!/bin/bash

#SBATCH --account=prj-smb018
#SBATCH --job-name=alternating-exponential
#SBATCH --output=alternating-exponential.out
#SBATCH --error=alternating-exponential.err
#SBATCH --time=02:00:00

julia --project=~/git/MassConversion.jl ~/git/MassConversion.jl/examples/alternating-exponential.jl
