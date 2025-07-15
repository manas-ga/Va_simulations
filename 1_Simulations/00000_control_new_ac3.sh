#!/bin/bash

#SBATCH --job-name=TEST
#SBATCH --output=/mnt/hel/obbard/Va_simulations/analyses/b_Interim_files/std_out/job_%A_%a.log   # Store logs in a custom directory
#SBATCH --open-mode=append                                                                       # Append output if the file already exists
#SBATCH --array=1-100%20                                                                         # Run replicate tasks
#SBATCH --ntasks=2
#SBATCH --mem=5G

# Path to Conda (if Conda isn't initialized by default)
CONDA_PATH="/home/msamant/miniconda3"  # Update to where Miniconda/Anaconda is installed
source $CONDA_PATH/etc/profile.d/conda.sh  # Initialize Conda in the job script

# Activate the Conda environment
conda activate marun

Param=`cat 000_parameter_grid_ac3.txt | awk "NR==${SLURM_ARRAY_TASK_ID}"`

Rscript 00_Control_sims.R $Param
