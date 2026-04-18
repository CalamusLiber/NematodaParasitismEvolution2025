#!/bin/sh

# ==============================================================================
# Trinity.DeConX.longestORF.sh
# -----------------------------------------------------------------------------
# Purpose: Extract the longest isoform per gene from decontaminated Trinity RNA-seq assemblies
#
# This script processes Trinity-assembled transcripts that have undergone decontamination
# (vector removal and taxonomic filtering). For each assembly, it selects the longest
# isoform sequence per Trinity gene using Trinity's built-in utility script.
#
# Input: Decontaminated FASTA files (*.DeContX.fasta) in the working directory
# Output: Longest isoform FASTA files (*_longest.fasta) in trans_ass_sift_longest/
# Logs: Progress and status messages written to ../DeContX.longestORF.log
#
# Dependencies: Trinity RNA-seq assembler, Perl, GNU parallel utilities
# Environment: HPC cluster with SLURM (24-thread parallel processing)
# ==============================================================================

# Load required HPC modules (some commented out as they may be pre-loaded or optional)
# module load lang/perl/5.30.0-bioperl-gcc 
# module load apps/trinity/2.8.5
module load apps/samtools/1.9
# module load apps/bowtie2/2.3.5
# module load lib/boost/1.71.0
module load lang/intel-parallel-studio-xe/2020

# Set up Python environment (using conda installation)
alias python='/To/Your/Directory/miniconda3/bin/python3.9'
alias python3='/To/Your/Directory/miniconda3/bin/python3.9'
export PATH=$PATH:/To/Your/Directory/miniconda3/bin:/To/Your/Directory/miniconda3/pkgs

# Set Trinity home directory for utility scripts
export TRINITY_HOME=/To/Your/Directory/trinityrnaseq

# Change to the directory containing decontaminated Trinity assemblies
cd $WORK/nematoda/data/rnaseq/trans_ass_sift

# Create output directory for longest isoform files
mkdir $WORK/nematoda/data/rnaseq/trans_ass_sift_longest

# Generate list of decontaminated assembly files (remove .fasta extension for processing)
rnaseq_list=$(ls *.DeContX* | sed -e "s/.fasta//g")

# ==============================================================================
# Multithread preparation using named pipe (FIFO) for parallel processing
# -----------------------------------------------------------------------------
# Create a temporary FIFO file for coordinating 24 parallel threads
# This allows controlled parallel execution without overloading the system
# ==============================================================================
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

# Initialize the FIFO with 24 tokens (one per thread)
for ((i=1;i<=24;i++))
do
    echo >&9
done

# ==============================================================================
# Process each decontaminated assembly in parallel
# -----------------------------------------------------------------------------
# For each assembly file, extract the longest isoform per Trinity gene
# Uses Trinity's get_longest_isoform_seq_per_trinity_gene.pl utility
# ==============================================================================
for nrna in $rnaseq_list
do
    # Acquire a thread token from the FIFO
    read -u9
    {
    
    # Check if the input FASTA file has content (size > 0 bytes)
    if [[ $(ls -l ./"$nrna".fasta | awk '{ print $5 }') -gt 0 ]]; then
        # Log the start of longest isoform selection
        echo -e "Selecting the longest isoform of "$nrna".fasta.\n" >> ../DeContX.longestORF.log
        
        # Execute Trinity's longest isoform extraction script
        # This Perl script parses Trinity gene/transcript headers and selects the longest sequence per gene
        perl ${TRINITY_HOME}/util/misc/get_longest_isoform_seq_per_trinity_gene.pl ./"$nrna".fasta > ../trans_ass_sift_longest/"$nrna"_longest.fasta && 
        
        # Log successful completion with output location
        echo -e "**********************************************\n" >> ../DeContX.longestORF.log
        echo -e ">>> Longest isoform of "$nrna".fasta is stored in trans_ass_sift_longest/"$nrna"_longest.fasta.\n" >> ../DeContX.longestORF.log
        echo -e "**********************************************\n" >> ../DeContX.longestORF.log
    else 
        # Log that the input file is empty and skip processing
        echo -e "**********************************************\n" >> ../DeContX.longestORF.log
        echo -e "But "$nrna".fasta file is empty.\n" >> ../DeContX.longestORF.log 
        echo -e "**********************************************\n\n" >> ../DeContX.longestORF.log
    fi
    
    # Release the thread token back to the FIFO
    echo >&9 
    } &   
done

# Wait for all background processes to complete
wait

