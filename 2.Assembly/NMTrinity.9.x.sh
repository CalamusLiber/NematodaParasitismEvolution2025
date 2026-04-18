#!/bin/sh
#SBATCH --job-name=NMTrinity.9

# ── Environment Setup ─────────────────────────────────────────────────────────

# Legacy module loads and aliases (currently disabled)
# module load lang/perl/5.30.0-bioperl-gcc 
# module load apps/trinity/2.8.5
# module load apps/samtools/1.9
# module load apps/bowtie2/2.3.5
# module load lib/boost/1.71.0
# module load lang/intel-parallel-studio-xe/2020
# alias python='/To/Your/Directory/anaconda3/bin/python3.9'
# alias python3='/To/Your/Directory/anaconda3/bin/python3.9'

# Update PATH to include Anaconda binaries and the specific RNA-seq environment
export PATH=$PATH:/To/Your/Directory/anaconda3/bin:/To/Your/Directory/anaconda3/envs/rnaseqenv/bin

# Legacy Trinity and TransDecoder environment variables (currently disabled)
# export TRINITY_HOME=/To/Your/Directory/trinityrnaseq
# export TDC_DIR=/To/Your/Directory/TransDecoder

# ── Directory Setup ───────────────────────────────────────────────────────────

# Navigate to the main RNA-seq working directory
cd /To/Your/Directory/nematoda/data/rnaseq

# Define base paths and working directories
DIR_cur=$(pwd)
HOME=/To/Your/Directory
DIR_cnd=/To/Your/Directory/anaconda3

# Activate the Conda environment for RNA-seq analysis
source $DIR_cnd/bin/activate $DIR_cnd/envs/rnaseqenv

# Define output directories for assembly and TransDecoder results
DIR_asm=$DIR_cur/trans_assembly
DIR_tdcout=$DIR_cur/transdecodes
DIR_tdcoutLg=$DIR_cur/transdecodesLongest

# Create necessary output directories if they do not exist
mkdir -p ./AfterQC.reports
mkdir -p $DIR_asm
mkdir -p $DIR_tdcout
mkdir -p $DIR_tdcoutLg

# ── Job Initialization ────────────────────────────────────────────────────────

# Get the sample name from the first command-line argument
nrna=$1

# Initialize the log file with job details and record the start time
echo -e "#######################\nTrinity Assembling start...\n$(date)\n#######################\n\n" >> NMTrinity.9.$nrna.log && 
echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> NMTrinity.9.$nrna.log && 
time0=$(date) && start_sec=$(date --date="$time0" +%s) 

# Define input FASTQ file paths
fq1=$DIR_cur/$nrna.R1.fastq
fq2=$DIR_cur/$nrna.R2.fastq

# ── Quality Control & Trimming ────────────────────────────────────────────────

# Run AfterQC for quality control and read trimming if the output directory doesn't exist
if [ ! -d $nrna.AfterQC.output ]; then 
    echo "Trimming $nrna starts." >> NMTrinity.9.$nrna.log

    mkdir $nrna.AfterQC.output	

    # Execute AfterQC using PyPy for performance, redirecting output to the log
    $HOME/pypy2/bin/pypy $HOME/AfterQC/after.py -1 $fq1 -2 $fq2 -g $nrna.AfterQC.output -b $nrna.AfterQC.output -r AfterQC.reports >> NMTrinity.9.$nrna.log 2>&1 && 

    # Remove bad quality reads to save disk space
    rm $nrna.AfterQC.output/*.bad.fq && 

    echo "Trimming $nrna is done." >> NMTrinity.9.$nrna.log
fi

# ── De Novo Assembly ──────────────────────────────────────────────────────────

# Run Trinity for de novo transcriptome assembly if the output FASTA doesn't exist
if [[ ! -f ./trans_assembly/"$nrna"_trinity.Trinity.fasta ]]; then
    # Execute Trinity with specified memory and CPU resources
    Trinity --seqType fq --max_memory 150G --left $nrna.AfterQC.output/$nrna.R1.good.fq --right $nrna.AfterQC.output/$nrna.R2.good.fq --CPU $[$SLURM_NTASKS-2] --output "$nrna"_trinity >> NMTrinity.9.$nrna.log 2>&1 &&
    
    # Move the resulting assembly files to the designated assembly directory
    mv "$nrna"_trinity.Trinity.fasta* ./trans_assembly && 
    echo "Assembling $nrna is done." >> NMTrinity.9.$nrna.log 
else
    echo ""$nrna"_trinity.Trinity.fasta exists, skip Trinity." >> NMTrinity.9.$nrna.log 
fi

# ── Post-Assembly Processing ──────────────────────────────────────────────────

# Check if the Trinity assembly was successfully generated and is not empty
if [[ $(ls -l ./trans_assembly/"$nrna"_trinity.Trinity.fasta | awk '{ print $5 }') -gt 0 ]]; then
    echo "Selecting the longest isoform." >> NMTrinity.9.$nrna.log && 
    
    # Extract the longest isoform per Trinity gene
    perl $DIR_cnd/envs/rnaseqenv/opt/trinity-2.13.2/util/misc/get_longest_isoform_seq_per_trinity_gene.pl ./trans_assembly/"$nrna"_trinity.Trinity.fasta > ./trans_assembly/"$nrna"_trinity_longest.fasta && 
    
    echo -e "**********************************************\nXXX>>> $nrna Trinity job is done successfully.\n**********************************************" >> NMTrinity.9.$nrna.log && 
    
    # Clean up intermediate AfterQC and Trinity directories to save space
    rm -R ./$nrna.AfterQC.output && 
    rm -R ./"$nrna"_trinity
else 
    echo -e "**********************************************\nXXXXX $nrna Trinity job expires.\nXXXXX But "$nrna"_trinity.Trinity.fasta file was not found.\n**********************************************" >> NMTrinity.9.$nrna.log
fi

# ── ORF Prediction (Full Assembly) ────────────────────────────────────────────

echo "TransDecoding $nrna starts." >> NMTrinity.9.$nrna.log

# Extract long open reading frames (ORFs) from the full Trinity assembly
TransDecoder.LongOrfs -t $DIR_asm/"$nrna"_trinity.Trinity.fasta --gene_trans_map $DIR_asm/"$nrna"_trinity.Trinity.fasta.gene_trans_map >> NMTrinity.9.$nrna.log 2>&1 &&

# Predict likely coding regions
TransDecoder.Predict -t $DIR_asm/"$nrna"_trinity.Trinity.fasta >> NMTrinity.9.$nrna.log 2>&1 &&

# Move TransDecoder output files to the designated directory
mv $DIR_cur/"$nrna"*.transdecoder.* $DIR_tdcout &&

# Rename the peptide sequence file for clarity
mv $DIR_tdcout/"$nrna"*.transdecoder.pep $DIR_tdcout/"$nrna"_trinity.tdc.pep.faa && 

# Clean up the TransDecoder temporary directory
rm -R $DIR_cur/$nrna*.transdecoder_dir* &&

# ── ORF Prediction (Longest Isoforms) ─────────────────────────────────────────

echo "TransDecoding $nrna LongestORF starts." >> NMTrinity.9.$nrna.log && 

# Extract long open reading frames (ORFs) from the longest isoforms
TransDecoder.LongOrfs -t $DIR_asm/"$nrna"_trinity_longest.fasta --gene_trans_map $DIR_asm/"$nrna"_trinity.Trinity.fasta.gene_trans_map >> NMTrinity.9.$nrna.log 2>&1 && 

# Predict likely coding regions for the longest isoforms
TransDecoder.Predict -t $DIR_asm/"$nrna"_trinity_longest.fasta >> NMTrinity.9.$nrna.log 2>&1 && 

# Move TransDecoder output files to the designated directory
mv $DIR_cur/"$nrna"*.transdecoder.* $DIR_tdcoutLg && 

# Rename the peptide sequence file for clarity
mv $DIR_tdcoutLg/"$nrna"*.transdecoder.pep $DIR_tdcoutLg/"$nrna"_trinity.tdclongest.pep.faa && 

# Clean up the TransDecoder temporary directory
rm -R $DIR_cur/$nrna*.transdecoder_dir* && 

echo "TransDecoding $nrna LongestORF finished." >> NMTrinity.9.$nrna.log 

# ── Job Completion & Timing ───────────────────────────────────────────────────

# Calculate total elapsed time and log the completion of the job
time1=$(date) && end_sec=$(date --date="$time1" +%s) && 
dsec=$((end_sec-start_sec)) && dhour=$(echo "scale=4; $dsec / 3600" | bc) && 
echo -e "#######################\nTransDecoding $nrna finished.\nAll assembling work accomplished.\n$time1\nOverall elapsed time: $dsec sec (= $dhour hr).\n#######################\n\n" >> NMTrinity.9.$nrna.log 
