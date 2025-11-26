#!/bin/bash
#SBATCH --job-name=mdp_ds0l
#SBATCH --account=xinyucl0
#SBATCH --mail-user=xinyucl@umich.edu
#SBATCH --time=168:00:00 # Request 168 hours of time
#SBATCH --mem-per-cpu=100000m # Request 100GB of memory per CPU
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

module load "python" # Load the correct module
# source ~/.bashrc && mamba activate <your_environment_name> # Activate custom environment if needed

# python runfile.py # The command to run your script
python valueIter.py