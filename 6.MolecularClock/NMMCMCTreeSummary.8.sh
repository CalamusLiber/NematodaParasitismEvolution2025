#!/bin/bash
# =============================================================================
# NMMCMCTreeSummary.8.sh - MCMCtree Results Summarization (Version 8)
# =============================================================================
#
# DESCRIPTION:
#   This script automates the summarization of completed MCMCtree analyses using
#   a SLURM-based workflow. It prepares and submits individual SLURM jobs for
#   each analysis combination to generate final reports and visualizations.
#
# PURPOSE:
#   Version 8 of the summarization script uses a more robust approach with
#   dedicated SLURM jobs for each analysis. This provides better resource
#   management, error handling, and parallel processing capabilities compared
#   to background execution in version 6.
#
# ANALYSIS GRID (Version 8 Configuration):
#   - Tree Calibrations: Automatically detected from .tre files in directory
#   - Clock Models: Independent Rates (IR) and Autocorrelated Rates (AR)
#   - Partitions: Single partition (1P) and five partitions (5P)
#   - Rate Heterogeneity: Gamma shape parameters (15, 20, 25)
#
# WORKFLOW:
#   1. Quality Control: Verify completion of all MCMC runs
#   2. Resource Allocation: Adjust memory based on partition complexity
#   3. Job Preparation: Copy scripts and configure SLURM parameters
#   4. Job Submission: Submit individual SLURM jobs for each analysis
#
# RESOURCE ALLOCATION:
#   - 5P analyses: 90GB RAM (more complex likelihood calculations)
#   - 1P analyses: 50GB RAM (simpler single partition)
#
# QUALITY CONTROL CRITERIA:
#   - No existing final report (avoids re-processing)
#   - Exactly 8 FigTree.tre files (4 posterior + 4 prior runs)
#   - Exactly 8 mcmc.txt files (4 posterior + 4 prior runs)
#
# OUTPUT FILES (per analysis, via SLURM job):
#   - MCMCtree_final_report.txt: Comprehensive analysis summary
#   - Time-calibrated phylogeny files
#   - MCMC diagnostics and trace plots
#   - Posterior vs prior distribution comparisons
#
# DEPENDENCIES:
#   - SumMCMCtree.py: Python script for MCMC summarization
#   - runPy.sh: SLURM wrapper script for Python execution
#   - Classification.xlsx: Taxon classification data
#   - tipcolors.xlsx: Color scheme for tree visualization
#   - Completed MCMCtree analyses with standard output files
#
# SLURM JOB CONFIGURATION:
#   - Automatic memory adjustment based on analysis complexity
#   - Individual jobs for each analysis combination
#   - Proper error handling and logging
# =============================================================================

# Get current working directory (root of analysis tree)
DIR_rt=$(pwd)

# Extract topology name from current directory path
Topo=$(echo $(pwd) | awk -v FS='/' '{print $NF}' | awk -v FS='.BD' '{print $1}')

# Auto-detect available tree calibration strategies from .tre files
TreeCalibration=$(ls *.tre | sed "s#\/##g;s#Nema.${Topo}.##g;s#.tre##g")

# Define analysis parameters for version 8 (full grid)
ClockModel='IR AR'        # Clock models: Independent and Autocorrelated Rates
Partition='1P 5P'         # Partition schemes: 1 and 5 partitions
rgene_gamma='15 20 25'    # Rate heterogeneity gamma shape parameters

# =============================================================================
# MAIN JOB SUBMISSION LOOP
# =============================================================================
# Prepare and submit SLURM jobs for each combination of analysis parameters
# to perform MCMCtree result summarization.

# Loop through tree calibration strategies
for TreCali in $TreeCalibration
do
    # Loop through partition schemes
    for Part in $Partition
    do
        # Loop through clock models
        for Clock in $ClockModel
        do
            # Loop through rate heterogeneity parameters
            for Rgene in $rgene_gamma
            do
                # Construct analysis directory path
                DIR_cur=$DIR_rt/${TreCali}.${Part}.${Clock}.rg${Rgene}

                # Navigate to analysis directory
                cd $DIR_cur

                # Check if summarization job needs to be submitted:
                # - No final report file exists (avoid re-processing)
                # - Exactly 8 FigTree.tre files (4 DIV + 4 PRD runs)
                # - Exactly 8 mcmc.txt files (4 DIV + 4 PRD runs)
                if [ $(find . -type f -name "*MCMCtree_final_report.txt" | wc -l) -lt 1 -a $(find . -type f -name FigTree.tre | wc -l) -eq 8 -a $(find . -type f -name mcmc.txt | wc -l) -eq 8 ]
                then
                    # Copy required scripts to analysis directory
                    cp /ampha/tenant/hebnu/private/user/hebnu204935/nematoda/timing/SumMCMCtree.py .
                    cp /ampha/tenant/hebnu/private/user/hebnu204935/nematoda/timing/runPy.sh .

                    # Make scripts executable
                    chmod -R 777 ./SumMCMCtree.py ./runPy.sh

                    # Generate output prefix from directory path
                    px=$(echo $DIR_cur | sed -s 's#/ampha/tenant/hebnu/private/user/hebnu204935/nematoda/timing/#NMtimetree_#g;s#\/#.#g')

                    # Adjust memory allocation based on partition complexity
                    # 5P analyses require more memory for complex likelihood calculations
                    if [ $Part == 1P ]
                    then
                        # Reduce memory for single partition analyses
                        sed -i 's/#SBATCH --mem-per-cpu=90GB/#SBATCH --mem-per-cpu=50GB/g' ./runPy.sh
                    fi

                    # Submit SLURM job for this analysis
                    # runPy.sh will execute SumMCMCtree.py with appropriate parameters
                    sbatch ./runPy.sh $DIR_cur $px
                fi
            done
        done
    done
done

# =============================================================================
# CLEANUP AND COMPLETION
# =============================================================================
# Return to original working directory and display completion message

# Return to the directory where the script was launched
cd $DIR_start

# Display completion message
echo "MCMCtree summarization jobs submitted successfully!"
echo "Each analysis combination will be processed in parallel via SLURM."
echo "Monitor job progress with 'squeue -u $USER' and check results in individual analysis directories."

