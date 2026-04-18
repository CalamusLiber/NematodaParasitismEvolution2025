#!/bin/bash
#SBATCH --job-name=NMTmTr2
#SBATCH --partition=HPC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB

# =============================================================================
# NMMCMCTr.8.x.2.sh - MCMCtree Divergence Time and Prior Estimation (Version 8)
# =============================================================================
#
# DESCRIPTION:
#   This SLURM job script performs the main MCMCtree analysis for molecular
#   clock dating, including divergence time estimation and prior distribution
#   sampling. It runs individual MCMC chains for Bayesian phylogenetic dating.
#   This is version 8 with enhanced quality control and job management.
#
# PURPOSE:
#   This script executes the core MCMC sampling phase of molecular clock dating
#   using MCMCtree. It performs two key analyses:
#   - Step 2: Estimates divergence times using sequence data
#   - Step 3: Samples from the prior distribution (for model comparison)
#
# PARAMETERS:
#   $1 (TreCali): Tree calibration strategy (e.g., CST1, CST2, etc.)
#   $2 (Part): Partition scheme (e.g., 1P, 5P for single/multiple partitions)
#   $3 (Clock): Clock model (IR=Independent Rates, AR=Autocorrelated Rates)
#   $4 (Rgene): Rate heterogeneity parameter for gamma distribution
#   $5 (R): Run number (for multiple independent MCMC chains)
#
# WORKFLOW:
#   1. MCMCtree Step 2: Estimate divergence times with sequence data
#   2. MCMCtree Step 3: Sample from prior distribution (no sequence data)
#   3. Enhanced quality control with line count verification
#   4. Automatic job resubmission for failed runs
#
# OUTPUT:
#   - DIR_DIV/run$R/: Divergence time estimation results
#   - DIR_PRD/run$R/: Prior distribution sampling results
#   - MCMC log files, tree files, and parameter estimates
#
# DEPENDENCIES:
#   - PAML software (mcmctree)
#   - Branch length estimates from step 1 (in.BV file)
#   - Control file template (mcmctree.ctl)
#
# SLURM RESOURCES:
#   - HPC partition, 1 node, 1 task, 1 CPU, 20GB RAM per CPU
#   - Enhanced job management with automatic resubmission
# =============================================================================

# Load PAML module and set executable paths
#module load paml/4.10.9
AppDir=/ampha/tenant/hebnu/private/user/hebnu204935/apps/paml-4.10.9/bin
MCMCTREE=$AppDir/mcmctree

# Command line arguments
TreCali=$1  # Tree calibration strategy
Part=$2     # Partition scheme
Clock=$3    # Clock model (IR/AR)
Rgene=$4    # Rate heterogeneity parameter
R=$5        # Run number for this MCMC chain
stopgen=3000000  # Target number of MCMC generations
sampfreq=15      # MCMC sampling frequency (every 15 generations)
DIR_rt=$(pwd)    # Root directory
DIR_cur=$DIR_rt/${TreCali}.${Part}.${Clock}.rg${Rgene}  # Analysis directory 

# =============================================================================
# MCMCtree Step 2: Estimate Divergence Times with Sequence Data
# =============================================================================
# This section performs the main Bayesian molecular clock analysis using
# sequence data to estimate divergence times. It uses the branch lengths
# estimated in Step 1 and applies fossil calibrations to date the phylogeny.
# Enhanced quality control includes verification of MCMC sample count.

# Navigate to analysis directory and create divergence time directory
cd $DIR_cur
mkdir -p $DIR_cur/DIR_DIV && cd $DIR_cur/DIR_DIV

# Check if Step 2 needs to be run (missing output, incomplete MCMC, or wrong sample count)
# Sample count check: Expected lines = (stopgen/sampfreq) + 2 (header + final line)
if [[ ! -f $DIR_cur/DIR_DIV/run$R/FigTree.tre ]] || [[ $(awk 'END{print $1}' $DIR_cur/DIR_DIV/run$R/mcmc.txt) -ne $stopgen ]] || [[ $(wc -l $DIR_cur/DIR_DIV/run$R/mcmc.txt | awk '{print $1}') -ne $[$stopgen/$sampfreq+2] ]]
then
    # Set up run directory
    mkdir -p $DIR_cur/DIR_DIV/run$R

    # Copy required files from branch length estimation
    cp $DIR_cur/DIR_BV/mcmctree.ctl $DIR_cur/DIR_DIV/run$R
    cp $DIR_cur/DIR_BV/in.BV $DIR_cur/DIR_DIV/run$R

    # Configure for divergence time estimation (usedata = 2)
    sed -i 's/usedata = 3/usedata = 2/g' $DIR_cur/DIR_DIV/run$R/mcmctree.ctl

    # Run MCMCtree Step 2
    cd $DIR_cur/DIR_DIV/run$R &&
    echo -e "#######################\nMCMCTREE STEP 2 start...\n$(date)\n#######################\n\n" > step2_log.text &&
    echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> step2_log.text &&
    mcmctree mcmctree.ctl >> step2_log.text 2>&1

    # Enhanced quality control: Resubmit job if analysis didn't complete successfully
    if [[ ! -f $DIR_cur/DIR_DIV/run$R/FigTree.tre ]] || [[ $(awk 'END{print $1}' $DIR_cur/DIR_DIV/run$R/mcmc.txt) -ne $stopgen ]] || [[ $(wc -l $DIR_cur/DIR_DIV/run$R/mcmc.txt | awk '{print $1}') -ne $[$stopgen/$sampfreq+2] ]]
    then
        rm -r $DIR_cur/DIR_DIV/run$R
        sbatch ./NMMCMCTr.8.x.2.sh $TreCali $Part $Clock $Rgene $R && scancel $SLURM_JOBID
    fi
fi 
		
# =============================================================================
# MCMCtree Step 3: Sample from Prior Distribution (No Sequence Data)
# =============================================================================
# This section samples from the prior distribution without using sequence data.
# This allows comparison between posterior and prior distributions to assess
# the information content of the molecular data and the influence of calibrations.
# Uses an until loop for automatic re-runs with enhanced quality control.

# Navigate to analysis directory and create prior distribution directory
cd $DIR_cur
mkdir -p $DIR_cur/DIR_PRD && cd $DIR_cur/DIR_PRD

# Check if Step 3 needs to be run (missing output, incomplete MCMC, wrong sample count, or seed mismatch)
if [[ ! -f $DIR_cur/DIR_PRD/run$R/FigTree.tre ]] || [[ $(awk 'END{print $1}' $DIR_cur/DIR_PRD/run$R/mcmc.txt) -ne $stopgen ]] || [[ $(wc -l $DIR_cur/DIR_PRD/run$R/mcmc.txt | awk '{print $1}') -ne $[$stopgen/$sampfreq+2] ]] || [[ $(cat $DIR_cur/DIR_DIV/run$R/SeedUsed) -ne $(grep -s 'seed =' $DIR_cur/DIR_PRD/run$R/mcmctree.ctl | awk -F ' = ' '{print $2}' | tr -d '\r') ]]
then
    # Clean up any existing run directory
    [[ -d $DIR_cur/DIR_PRD/run$R ]] && rm -r $DIR_cur/DIR_PRD/run$R

    # Set up fresh run directory
    mkdir -p $DIR_cur/DIR_PRD/run$R

    # Copy control file template
    cp $DIR_cur/DIR_BV/mcmctree.ctl $DIR_cur/DIR_PRD/run$R

    # Configure for prior sampling (usedata = 0, no sequence data)
    sed -i 's/usedata = 3/usedata = 0/g' $DIR_cur/DIR_PRD/run$R/mcmctree.ctl

    # Change output file name to distinguish from posterior
    sed -i 's/outfile = out/outfile = out_prd/g' $DIR_cur/DIR_PRD/run$R/mcmctree.ctl

    # Use same random seed as Step 2 for reproducible comparison
    sed -i "s#seed = -1#seed = $(cat $DIR_cur/DIR_DIV/run$R/SeedUsed)#g" $DIR_cur/DIR_PRD/run$R/mcmctree.ctl

    # Run MCMCtree Step 3
    cd $DIR_cur/DIR_PRD/run$R &&
    echo -e "#######################\nMCMCTREE STEP 3 start...\n$(date)\n#######################\n\n" > step3_log.text &&
    echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> step3_log.text &&
    mcmctree mcmctree.ctl >> step3_log.text 2>&1

    # Enhanced quality control: Automatic re-run until convergence criteria met
    until [[ -f $DIR_cur/DIR_PRD/run$R/FigTree.tre ]] && [[ $(awk 'END{print $1}' $DIR_cur/DIR_PRD/run$R/mcmc.txt) -ge $stopgen ]] && [[ $(wc -l $DIR_cur/DIR_PRD/run$R/mcmc.txt | awk '{print $1}') -eq $[$stopgen/$sampfreq+2] ]]
    do
        # Clean up incomplete MCMC output and re-run
        rm $DIR_cur/DIR_PRD/run$R/mcmc.txt
        echo -e "#######################\nMCMCTREE STEP 3 start...\n$(date)\n#######################\n\n" > step3_log.text &&
        echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> step3_log.text &&
        mcmctree mcmctree.ctl >> step3_log.text 2>&1
    done
fi


# =============================================================================
# MCMCtree Step 4: Summarize MCMC Results (Deprecated/Commented Out)
# =============================================================================
# This section would summarize the MCMC results from multiple runs, combining
# posterior and prior distributions for final analysis and visualization.
# Currently commented out - would require SumMCMCtree.py script and data files.
#
# Functionality would include:
# - Combine results from multiple MCMC chains
# - Generate summary statistics and plots
# - Create time-calibrated phylogenies
# - Compare posterior vs prior distributions
#
# Note: This step is typically run separately after all MCMC chains complete.

# cd $DIR_cur
# cp /To/Your/Directory/nematoda/SumMCMCtree.py .
# chmod -R 777 ./SumMCMCtree.py
# px=$(echo $DIR_cur | sed -s 's#/To/Your/Directory/nematoda/#NMtimetree_#g;s#\/#.#g')
# ./SumMCMCtree.py --postdir $DIR_cur/DIR_DIV --priodir $DIR_cur/DIR_PRD --classification /To/Your/Directory/nematoda/Classification.xlsx --tipcolor /To/Your/Directory/nematoda/tipcolors.xlsx --prefix $px 
