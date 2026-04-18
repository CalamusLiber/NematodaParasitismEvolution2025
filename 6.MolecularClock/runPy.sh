#!/bin/bash
#SBATCH --job-name=SumMCMCTr
#SBATCH --partition=HPC-HT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90GB
#SBATCH --time=5:00:0

# runPy.sh - SLURM wrapper script for MCMCtree result summarization
# This script executes SumMCMCtree.py on HPC systems where direct Python execution
# via SLURM is not supported. It handles environment setup and parameter passing.
#
# Usage: sbatch runPy.sh <analysis_directory> <output_prefix>
# Parameters:
#   $1 (DIR_cur): Path to the analysis directory containing DIR_DIV and DIR_PRD
#   $2 (px): Output prefix for generated files
#
# Dependencies:
#   - SumMCMCtree.py: Python script for MCMCtree output processing
#   - Classification.xlsx: Taxonomic classification data
#   - tipcolors.xlsx: Color scheme for tree tips
#   - Anaconda3 module with required Python packages

# Specify Python executable path for systems without default Python 3.12
PYTHON=/To/your/python/directory/python3.12

# Command line arguments
DIR_cur=$1  # Analysis directory containing MCMC output folders
px=$2       # Output file prefix

# Navigate to the analysis directory
cd $DIR_cur

# Load required modules (Anaconda3 for Python environment)
module load anaconda3/2024.06

# Execute SumMCMCtree.py with specified parameters:
# --postdir: Directory containing posterior MCMC samples (DIR_DIV)
# --priodir: Directory containing prior MCMC samples (DIR_PRD)
# --classification: Excel file with taxonomic classifications for coloring
# --tipcolor: Excel file with color assignments for tree tips
# --prefix: Output file prefix
# --burnin: Number of samples to discard as burn-in (50000)
# --select_chains: Automatically select best MCMC chains
# --good_ess: Minimum effective sample size threshold (100)
# --rm_mcmc: MCMC file handling ('none' = keep original files)
# --compress: Compress output files
$PYTHON ./SumMCMCtree.py --postdir $DIR_cur/DIR_DIV --priodir $DIR_cur/DIR_PRD --classification /ampha/tenant/hebnu/private/user/hebnu204935/nematoda/timing/Classification.xlsx --tipcolor /ampha/tenant/hebnu/private/user/hebnu204935/nematoda/timing/tipcolors.xlsx --prefix $px --burnin 50000 --select_chains --good_ess 100 --rm_mcmc none --compress 
