#!/bin/bash
# =============================================================================
# NMMCMCTreeChecking.sh - MCMCtree Analysis Quality Control and Status Check
# =============================================================================
#
# DESCRIPTION:
#   This script performs quality control checks on completed MCMCtree analyses.
#   It scans through all analysis directories to identify successfully completed
#   runs and detect potential issues with MCMC sampling consistency.
#
# PURPOSE:
#   After running multiple MCMCtree analyses with different calibration strategies,
#   clock models, partitions, and rate heterogeneity parameters, this script
#   helps identify which analyses completed successfully and flags any issues
#   that may require re-running or further investigation.
#
# ANALYSIS GRID:
#   - Tree Calibrations: Automatically detected from .tre files in directory
#   - Clock Models: Independent Rates (IR) and Autocorrelated Rates (AR)
#   - Partitions: Single partition (1P) and five partitions (5P)
#   - Rate Heterogeneity: Gamma shape parameters (15, 20, 25)
#   - MCMC Runs: 4 independent chains per analysis
#
# QUALITY CONTROL CHECKS:
#   1. Completion Status: Verifies presence of required output files
#   2. Sample Count Consistency: Compares MCMC sample counts between posterior
#      and prior runs to detect sampling issues
#   3. File Count Verification: Ensures all expected output files are present
#
# OUTPUT:
#   - Lists directories with completed analyses
#   - Reports any inconsistencies in MCMC sample counts between DIV and PRD runs
#   - Silent for analyses that pass all checks
#
# USAGE:
#   Run this script from the root directory containing all MCMCtree analysis
#   subdirectories. It will automatically scan and report on analysis status.
#
# DEPENDENCIES:
#   - Completed MCMCtree analyses with standard directory structure
#   - Standard output files: FigTree.tre, mcmc.txt, MCMCtree_final_report.txt
# =============================================================================

# Get current working directory (root of analysis tree)
DIR_rt=$(pwd)

# Extract topology name from current directory path
Topo=$(echo $(pwd) | awk -v FS='/' '{print $NF}' | awk -v FS='.BD' '{print $1}')

# Auto-detect available tree calibration strategies from .tre files
TreeCalibration=$(ls *.tre | sed "s#\/##g;s#Nema.${Topo}.##g;s#.tre##g")

# Define analysis parameters to check
ClockModel='IR AR'          # Clock models: Independent and Autocorrelated Rates
Partition='1P 5P'           # Partition schemes: 1 and 5 partitions
rgene_gamma='15 20 25'      # Rate heterogeneity gamma shape parameters
nruns=4                     # Number of independent MCMC runs per analysis

# =============================================================================
# MAIN ANALYSIS CHECKING LOOP
# =============================================================================
# Iterate through all combinations of analysis parameters to check completion status
# and quality control for each MCMCtree analysis.

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

                # Check completion criteria:
                # - No final report file (indicates summarization not done)
                # - Exactly 8 FigTree.tre files (4 DIV + 4 PRD runs)
                # - Exactly 8 mcmc.txt files (4 DIV + 4 PRD runs)
                if [ $(find . -type f -name "*MCMCtree_final_report.txt" | wc -l) -lt 1 -a $(find . -type f -name FigTree.tre | wc -l) -eq 8 -a $(find . -type f -name mcmc.txt | wc -l) -eq 8 ]
                then
                    # Report directory with completed MCMC runs
                    echo $DIR_cur

                    # Check for MCMC sampling consistency between posterior and prior runs
                    for ((R=1;R<=$nruns;R++))
                    do
                        # Count lines in posterior MCMC output (DIV)
                        div=$(wc -l ./DIR_DIV/run$R/mcmc.txt | awk '{print $1}')

                        # Count lines in prior MCMC output (PRD)
                        prd=$(wc -l ./DIR_PRD/run$R/mcmc.txt | awk '{print $1}')

                        # Report inconsistency if sample counts don't match
                        [ $div -ne $prd ] && echo -e "$div $DIR_cur/DIR_DIV/run$R\n$prd $DIR_cur/DIR_PRD/run$R"
                    done
                fi
            done
        done
    done
done

