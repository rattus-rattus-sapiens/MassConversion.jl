#!/bin/bash

#SBATCH --account=prj-smb018
#SBATCH --job-name=runsim
#SBATCH --output=runsim.out
#SBATCH --error=runsim.err
#SBATCH --time=02:00:00

julia --project=~/git/MassConversion.jl ~/git/MassConversion.jl/examples/runsim.jl