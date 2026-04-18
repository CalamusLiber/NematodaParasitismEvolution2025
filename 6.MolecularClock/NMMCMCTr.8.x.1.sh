#!/bin/bash
#SBATCH --job-name=NMTmTr1
#SBATCH --partition=HPC-HT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=30GB

# =============================================================================
# NMMCMCTr.8.x.1.sh - MCMCtree Analysis Setup and Initial Run (Version 8)
# =============================================================================
#
# DESCRIPTION:
#   This SLURM job script sets up and initiates MCMCtree analysis for molecular
#   clock dating of nematode phylogenies. It performs the initial setup phase
#   including branch length estimation and prepares for multiple MCMC runs.
#   This is version 8 with updated resource allocations and module loading.
#
# PURPOSE:
#   MCMCtree (from PAML package) implements Bayesian molecular clock dating.
#   This script handles the preparatory steps before running multiple MCMC
#   chains for robust phylogenetic dating with fossil calibrations.
#
# PARAMETERS:
#   $1 (TreCali): Tree calibration strategy (e.g., CST1, CST2, etc.)
#   $2 (Part): Partition scheme (e.g., 1P, 5P for single/multiple partitions)
#   $3 (Clock): Clock model (IR=Independent Rates, AR=Autocorrelated Rates)
#   $4 (Rgene): Rate heterogeneity parameter for gamma distribution
#
# WORKFLOW:
#   1. Set up working directory structure
#   2. Generate branch length estimates (in.BV file)
#   3. Launch multiple MCMC runs via SLURM
#
# OUTPUT:
#   - DIR_BV/: Contains branch length estimation results
#   - DIR_DIV/ and DIR_PRD/: Directories for MCMC runs (created by step 2)
#   - Log files and intermediate results
#
# DEPENDENCIES:
#   - PAML software (mcmctree, codeml, baseml) via module or AppDir
#   - Input files: sequence alignments, tree topologies, calibration priors
#   - mcmctree.ctl control file template
#
# SLURM RESOURCES:
#   - HPC-HT partition, 1 node, 1 task, 2 CPUs, 30GB RAM per CPU
#   - Uses sbatch for job submission instead of background processes
# =============================================================================

# Load PAML module and set executable paths
#module load paml/4.10.9
AppDir=/ampha/tenant/hebnu/private/user/hebnu204935/apps/paml-4.10.9/bin
MCMCTREE=$AppDir/mcmctree
BASEML=$AppDir/baseml
CODEML=$AppDir/codeml

# Get current working directory and extract topology name from path
DIR_rt=$(pwd)
cd $DIR_rt
Topo=$(echo $(pwd) | awk -v FS='/' '{print $NF}' | awk -v FS='.BD' '{print $1}')

# Command line arguments
TreCali=$1  # Tree calibration strategy
Part=$2     # Partition scheme
Clock=$3    # Clock model (IR/AR)
Rgene=$4    # Rate heterogeneity parameter
stopgen=3000000  # Target number of MCMC generations

# Create analysis directory structure
mkdir -p $DIR_rt/${TreCali}.${Part}.${Clock}.rg${Rgene} &&
cd $DIR_rt/${TreCali}.${Part}.${Clock}.rg${Rgene}
DIR_cur=$(pwd)
nruns=4  # Number of independent MCMC runs for convergence assessment

# =============================================================================
# MCMCtree Step 1: Generate Branch Length Estimates (in.BV file)
# =============================================================================
# This section estimates branch lengths using maximum likelihood under the
# specified substitution model. The in.BV file contains the estimated branch
# lengths that will be used as priors in the subsequent MCMC analysis.

# Check if branch length estimation is needed
if [[ ! -f $DIR_cur/DIR_BV/in.BV ]]
then
    cd $DIR_cur &&

    # Create directory for branch length estimation
    mkdir -p $DIR_cur/DIR_BV
    cd $DIR_cur/DIR_BV

    # Copy control file template
    cp $DIR_rt/mcmctree.ctl ./

    # Configure sequence file path
    sed -i "s#seqfile = seq.phy#seqfile = $DIR_rt/Nema.TBDL.nom.rmsp.${Part}.phy#" mcmctree.ctl

    # Configure tree file path (uses Topo variable for topology name)
    sed -i "s#treefile = input.tre#treefile = $DIR_rt/Nema.${Topo}.${TreCali}.tre#g" mcmctree.ctl

    # Optional: Set root age for specific calibration strategies
    # if [[ $TreCali == 'CST3' ]]
    # then
    #     sed -i "s#RootAge =#RootAge = 6.09#g" mcmctree.ctl
    # fi

    # Configure number of data partitions
    if [[ $Part == '5P' ]]
    then
        sed -i "s#ndata = 1#ndata = 5#g" mcmctree.ctl
    fi

    # Configure clock model
    if [[ $Clock == 'IR' ]]
    then
        sed -i "s#clock = #&2#g" mcmctree.ctl  # Independent rates
    elif [[ $Clock == 'AR' ]]
    then
        sed -i "s#clock = #&3#g" mcmctree.ctl  # Autocorrelated rates
    fi

    # Configure rate heterogeneity among sites
    sed -i "s#rgene_gamma = 2 x#rgene_gamma = 2 ${Rgene}#g" mcmctree.ctl

    # Start logging
    echo -e "#######################\nMCMCTREE STEP 1 start...\n$(date)\n#######################\n\n" > step1_log.text &&
    echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> step1_log.text &&

    # Run MCMCtree for initial branch length estimation
    mcmctree mcmctree.ctl >> step1_log.text 2>&1 &&

    echo -e "\n-----------mcmctree job done!-----------\n" >> step1_log.text &&

    # Copy amino acid rate matrix file
    cp $DIR_rt/lg.dat ./

    # Configure CODEML control files for each partition
    # Set model to empirical amino acid model, add rate matrix, and configure gamma
    for file in tmp*ctl
    do
        sed -i "s/model = 0/model = 2/g" $file           # Empirical AA model
        sed -i "s/aaRatefile = /aaRatefile = lg.dat/g" $file  # LG matrix
        sed -i "s/method = 0/method = 1/g" $file         # ML estimation
        echo -e "fix_alpha = 0\nalpha = 0.5\nncatG = 4" >> $file  # Gamma rates
    done

    # Clean up any existing in.BV file
    #[[ -f in.BV ]] && rm in.BV

    # Run CODEML for each partition to estimate branch lengths
    for file in tmp*ctl
    do
        codeml $file >> step1_log.text 2>&1 &&

        # Optional: Test for large branch lengths (commented out)
        # while [[ $(sed -n '6p' rst2 | awk '{for (i=1; i<=NF; i++) if ($i >= 20) {print "True"; exit}}') == 'True' ]]
        # do
        #     codeml $file
        # done

        # Clean up spacing in branch length output to prevent parsing errors
        sed -i '6s/[[:space:]]\{1,\}/  /g' rst2 &&
        sed -i '6s/[[:space:]]\{3,\}/  /g' rst2 &&

        # Append results to in.BV file
        cat rst2 >> in.BV
    done

    echo -e "\n-----------codeml job done!-----------\n" >> step1_log.text
fi 

# =============================================================================
# MCMCtree Step 2: Launch Multiple MCMC Runs via SLURM
# =============================================================================
# This section launches multiple independent MCMC chains for Bayesian molecular
# clock dating. Running multiple chains allows assessment of convergence and
# provides more robust posterior estimates. Uses sbatch for proper job queuing.

# Launch multiple MCMC runs via SLURM job submission
for ((R=1;R<=$nruns;R++))
do
    cd $DIR_rt

    # Check if MCMC run needs to be (re)started
    # Conditions: Missing output files OR incomplete runs (not exactly stopgen generations)
    if [[ ! -f $DIR_cur/DIR_DIV/run$R/FigTree.tre ]] || [[ ! -f $DIR_cur/DIR_PRD/run$R/FigTree.tre ]] || [[ $(awk 'END{print $1}' $DIR_cur/DIR_DIV/run$R/mcmc.txt) -ne $stopgen ]] || [[ $(awk 'END{print $1}' $DIR_cur/DIR_PRD/run$R/mcmc.txt) -ne $stopgen ]]
    then
        # Submit MCMC run as separate SLURM job
        sbatch ./NMMCMCTr.8.x.2.sh $TreCali $Part $Clock $Rgene $R
    fi
done

