#!/bin/bash
#SBATCH --job-name=NMBU.x          # Job name identifier
#SBATCH --exclusive                # Exclusive node access (no sharing)
#SBATCH --partition=mwvdk          # Partition/queue to submit to
#SBATCH --nodes=4                  # Number of compute nodes
#SBATCH --ntasks-per-node=4        # Tasks per node
#SBATCH --cpus-per-task=8          # CPU cores per task
#SBATCH --mem=360GB                # Total memory allocation

# ==============================================================================
# NMAlignTrim.sh - BUSCO Multiple Sequence Alignment, Trimming, and Quality Assessment
# -----------------------------------------------------------------------------
# Purpose: Process BUSCO gene sequences through alignment, trimming, and quality evaluation
#
# This SLURM script performs parallel processing of BUSCO genes with the following workflow:
# 1. Multiple sequence alignment using MAFFT (E-INS-i algorithm)
# 2. Alignment trimming using BMGE with two strategies:
#    - Normal trimming (stationary composition assumed)
#    - Stationary-based trimming (accounts for compositional heterogeneity)
# 3. Quality assessment using PhyKit to calculate various metrics
#
# Input: Directory with BUSCO FASTA files (*.faa) from UniteBUSCOs.py
# Output: Aligned and trimmed sequences, quality metrics CSV files
#
# Tools: MAFFT, BMGE, PhyKit
# Environment: HPC cluster with SLURM job scheduler
# ==============================================================================

# Calculate parallel processing parameters from SLURM environment
NP=$SLURM_CPUS_PER_TASK           # CPUs per task (8)
NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE]  # Total tasks across nodes (16)

# Load required software modules
# module load apps/iqtree/2.1.3    # IQ-TREE (commented out, not used in this script)
module load apps/bmge/1.12/       # BMGE alignment trimming tool
# module load apps/bmge/1.12      # Alternative BMGE version (commented out)

# Define tool executable paths
IQTREE=/sw/apps/iqtree-2.1.3-Linux/bin/iqtree2  # IQ-TREE binary (not used)
BMGE=/sw/apps/BMGE-1.12/BMGE.jar               # BMGE JAR file
# FASconCATg=/To/Your/Directory/FASconCAT-G/FASconCAT-G_v1.05.pl  # Concatenation tool (commented out)
# Define working directories
DIR_cnd=/To/Your/Directory/anaconda3                    # Conda installation directory
DIR_united=/To/Your/Directory/nematoda/data/busco/united  # Input: BUSCO united sequences
DIR_msa=/To/Your/Directory/nematoda/msa                 # Output: MSA processing directory
DIR_align=$DIR_msa/align                                # Output: Aligned sequences
DIR_trim=$DIR_msa/trim                                  # Output: Trimmed sequences
DIR_matrix=$DIR_msa/matrix                              # Output: Quality metrics
DIR_log=$DIR_msa/log                                    # Output: Processing logs

# Create output directory structure
mkdir -p $DIR_msa
mkdir -p $DIR_align
mkdir -p $DIR_trim
mkdir -p $DIR_trim/sbt    # Stationary-based trimming results
mkdir -p $DIR_trim/nom    # Normal trimming results
mkdir -p $DIR_matrix
mkdir -p $DIR_log

# Activate conda environment with required tools (MAFFT, PhyKit, etc.)
source $DIR_cnd/bin/activate $DIR_cnd/envs/msaenv

# ==============================================================================
# Stage 1: Multiple Sequence Alignment with MAFFT
# -----------------------------------------------------------------------------
# Align BUSCO protein sequences using MAFFT's E-INS-i algorithm
# This is suitable for sequences with global homology (typical for BUSCOs)
# ==============================================================================

cd $DIR_united && 
UniteList=$(ls *.faa | sed 's/.faa//g')  # Get list of BUSCO genes to process

# Set up parallel processing with named pipe (FIFO) for 16 concurrent tasks
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2  # Clean up FIFO on script exit
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

# Initialize FIFO with 16 tokens (one per parallel task)
for ((i=1;i<=$NT;i++))
do
    echo >&9
done

for gene in $UniteList
do 
    read -u9
    {

    if [[ ! -f $DIR_align/$gene.mft.faa ]] || [[ $(ls -l $DIR_align/$gene.mft.faa | awk '{print $5}') -eq 0 ]]; then 
        echo -e "#######################\nAligning, Trimming, and Filtering - 1\n$(date)\n#######################\n\n" >> $DIR_log/$gene.msatreat.log && 
        echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> $DIR_log/$gene.msatreat.log && 
        echo -e "=======================\nAligning with MAFFT start...\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log && 
        einsi --maxiterate 1000 --reorder --anysymbol --thread $NP $DIR_united/$gene.faa > $DIR_align/$gene.mft.faa && 
        echo -e "=======================\nAligning with MAFFT accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log 
    fi

    echo >&9 
    } & 
done
wait 

# ==============================================================================
# Stage 2: Alignment Trimming and Quality Assessment
# -----------------------------------------------------------------------------
# Trim poorly aligned regions using BMGE and assess alignment quality with PhyKit
# Two trimming strategies: normal (assumes stationary composition) and stationary-based
# ==============================================================================

cd $DIR_trim

tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile
for ((i=1;i<=$NT;i++))
do
    echo >&9
done

for gene in $UniteList
do 
    # Acquire processing token
    read -u9
    {
    
    # ==========================================================================
    # Normal BMGE trimming (assumes stationary amino acid composition)
    # ==========================================================================
    echo -e "=======================\nTrimming with BMGE (normal) start...\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log && 
    java -Xmx90G -jar $BMGE -i $DIR_align/$gene.mft.faa -t AA -m BLOSUM90 -h 0.4 -s NO -of $DIR_trim/nom/$gene.mft.bmge.nom.faa >> $DIR_log/$gene.msatreat.log 2>&1 && 
    echo -e "=======================\nTrimming with BMGE (normal) accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log && 

    # Assess quality metrics for normal-trimmed alignment
    echo -e "=======================\nFiltering with PhyKit (for normal BMGE output) start...\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log && 
    ptaxa1=$(grep '^>' $DIR_trim/nom/$gene.mft.bmge.nom.faa | wc -l) &&                                    # Count taxa/sequences
    parsisite1=($(phykit parsimony_informative_sites $DIR_trim/nom/$gene.mft.bmge.nom.faa)) &&           # Parsimony-informative sites (returns 3 values)
    GCcont1=$(phykit gc_content $DIR_trim/nom/$gene.mft.bmge.nom.faa) &&                                 # GC content
    RCV1=$(phykit relative_composition_variability $DIR_trim/nom/$gene.mft.bmge.nom.faa) &&             # Relative composition variability
    EvolRate1=$(phykit pairwise_identity $DIR_trim/nom/$gene.mft.bmge.nom.faa | sed -n 1p | awk '{print $2}') &&  # Evolutionary rate (from pairwise identity)
    
    # Write quality metrics to CSV file
    echo -e "$DIR_trim/nom/$gene.mft.bmge.nom.faa,$ptaxa1,${parsisite1[0]},${parsisite1[1]},${parsisite1[2]},$GCcont1,$RCV1,$EvolRate1" >> $DIR_matrix/LociFilterMeasure.nom.list && 

    echo -e "=======================\nFiltering with PhyKit (for normal BMGE output) accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log    
        
    echo >&9 
    } & 
done
wait 

tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile
for ((i=1;i<=$NT;i++))
do
    echo >&9
done

for gene in $UniteList
do 
    read -u9
    {
    
    echo -e "=======================\nTrimming with BMGE (stationary-based) start...\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log && 
    java -Xmx90G -jar $BMGE -i $DIR_align/$gene.mft.faa -t AA -m BLOSUM90 -h 0.4 -s YES -of $DIR_trim/sbt/$gene.mft.bmge.sbt.faa >> $DIR_log/$gene.msatreat.log 2>&1 && 
    echo -e "=======================\nTrimming with BMGE (stationary-based) accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log && 
    
    echo -e "=======================\nFiltering with PhyKit (for stationary-based BMGE output) start...\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log && 
    ptaxa2=$(grep '^>' $DIR_trim/sbt/$gene.mft.bmge.sbt.faa | wc -l) && 
    parsisite2=($(phykit parsimony_informative_sites $DIR_trim/sbt/$gene.mft.bmge.sbt.faa)) && # 3 values
    GCcont2=$(phykit gc_content $DIR_trim/sbt/$gene.mft.bmge.sbt.faa) && 
    RCV2=$(phykit relative_composition_variability $DIR_trim/sbt/$gene.mft.bmge.sbt.faa) && 
    EvolRate2=$(phykit pairwise_identity $DIR_trim/sbt/$gene.mft.bmge.sbt.faa | sed -n 1p | awk '{print $2}') && 
    echo -e "$DIR_trim/sbt/$gene.mft.bmge.sbt.faa,$ptaxa2,${parsisite2[0]},${parsisite2[1]},${parsisite2[2]},$GCcont2,$RCV2,$EvolRate2" >> $DIR_matrix/LociFilterMeasure.sbt.list && 

    echo -e "=======================\nFiltering with PhyKit (for stationary-based BMGE output) accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/$gene.msatreat.log &&         
    
    time1=$(date) && end_sec=$(date --date="$time1" +%s) && 
    dsec=$((end_sec-start_sec)) && dhour=$(echo "scale=4; $dsec / 3600" | bc) 

    echo -e "#######################\nAligning, Trimming, and Filtering - 1 accomplished.\n$time1\nOverall elapsed time: $dsec sec (= $dhour hr).\n#######################\n\n" >> $DIR_log/$gene.msatreat.log

    # Release processing token
    echo >&9 
    } & 
done
wait  # Wait for all processing to complete

