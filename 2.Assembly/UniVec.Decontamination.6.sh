#!/bin/sh
#SBATCH --job-name=NMDecont.5
#SBATCH --partition=compute
#SBATCH --nodes=4                   # Number of compute nodes
#SBATCH --ntasks-per-node=3         # Parallel jobs per node
#SBATCH --cpus-per-task=8           # CPU cores per job
#SBATCH --mem=350GB                 # Total memory allocation
#SBATCH --exclusive                 # Exclusive node access

# ============================================================================
# UniVec and Taxonomic Decontamination Pipeline for RNA-seq Assemblies
# ============================================================================
#
# This script performs a multi-stage decontamination workflow:
# 1. Vector removal: Removes cloning vector sequences using blastn against UniVec
# 2. Alien sequence detection: Uses DIAMOND blastx to identify non-target taxa
# 3. Taxonomic filtering: Removes sequences not matching the target organism's phylum
# 4. ORF prediction: Translates cleaned sequences using TransDecoder
# 5. Visualization: Generates pie charts showing transcript taxonomic composition

# Calculate parallelization parameters based on SLURM allocation
DNP=$SLURM_CPUS_PER_TASK
NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE]
FNP=$[$SLURM_NTASKS-2]

cd /To/Your/Directory/nematoda/data/rnaseq
DIR_cur=$(pwd)
DIR_cnd=/To/Your/Directory/anaconda3                   # Conda environment root
UniVecDB=/To/Your/Directory/db/UniVec/UniVec          # BLAST database of cloning vectors
NRDiaDB=/To/Your/Directory/db/nr/nr.dmnd              # DIAMOND database of NR proteins
SeqOrinPiePlot=/To/Your/Directory/nematoda/data/SeqOriginPiePlot.py  # Visualization script
TaxDumpDB=/To/Your/Directory/db/taxdump               # Local NCBI taxonomy dump (use "remote" for online queries, not recommended on HPC)

# Output directory structure
DIR_asm=$DIR_cur/trans_assembly                        # Input assemblies
DIR_sif=$DIR_cur/trans_ass_sift                        # Sifted (vector-removed) sequences
DIR_dec=$DIR_cur/decontamination                       # Decontamination workflow outputs
DIR_tdc=$DIR_cur/transdecodesX                         # TransDecoder outputs
DIR_Trc=$DIR_cur/TranscriptContents                    # Taxonomy composition visualizations

# Create output directories
mkdir -p $DIR_sif
mkdir -p $DIR_dec
mkdir -p $DIR_tdc
mkdir -p $DIR_Trc

# Identify input files and classify by source organism
InputFasta=$(ls ./trans_assembly/*_trinity.Trinity.fasta | sed -e "s/.fasta//g;s/.\/trans_assembly\///g")  # All input FASTA files
NNN=$(echo $InputFasta | wc -w)                        # Count of input files

# Identify Ecdysozoans that require special decontamination logic
# (Onychophora, Kinorhyncha, Tardigrada, Priapulida, Loricifera) validated against Ecdysozoan genes
EcdysFasta=$(ls ./trans_assembly/*_trinity.Trinity.fasta | sed -e "s/.fasta//g;s/.\/trans_assembly\///g" | grep -E "^OP_|^KR_|^TG_|^PA_|^LR_")

# ============================================================================
# Function Definitions
# ============================================================================

# runDeVect: Iteratively remove vector sequences using blastn and DeContaminRNA
# Performs multiple rounds of screening until no vectors are detected.
# Input: $1 = full path to FASTA file
runDeVect(){

    local rnatri=$(echo $1 | sed -e "s#$DIR_asm\/##g;s#.fasta##g")

    # Select decontamination program based on organism type
    # Ecdysozoans use DeContaminRNA-tx.py (taxonomy-aware), others use DeContaminRNA-rk.py (vector-focused)
    if [[ $EcdysFasta =~ $rnatri ]]; then 
        local DeConProg=/To/Your/Directory/nematoda/data/DeContaminRNA-tx.py
    else
        local DeConProg=/To/Your/Directory/nematoda/data/DeContaminRNA-rk.py
    fi
    
    echo -e "\n########################\nRemoval of vector sequences...\n$(date)\n########################\n" 
    echo -e "\n========================\nScreening $rnatri.fasta (round 1)...\n========================\n" 

    # Initialize with input assembly
    cp $1 $DIR_sif/$rnatri.Sifted.fasta
        
    echo -e "\nStart screening $rnatri.fasta.\n" 
    echo -e "\n$(blastn -version)\n" 
    # First BLASTN round against UniVec with optimized parameters for vector detection
    blastn -db $UniVecDB -query $DIR_sif/$rnatri.Sifted.fasta -out $DIR_dec/$rnatri.decont -outfmt "6 qseqid sseqid qstart qend staxids" -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 && 
    echo -e "\nFinished screening $rnatri.fasta.\n$(date)\n" 
    
    # Iteratively remove vectors and re-screen until no hits or max rounds reached
    local counter=1

    while [ $(ls -l $DIR_dec/$rnatri.decont | awk '{ print $5 }') -gt 0 -a $counter -le 10 ] 
    do
        echo -e "\n========================\nScreening $rnatri.fasta (round $[$counter+1])...\n$(date)\n========================\n" 
        echo -e "\nRemoving vector sequences from $rnatri.fasta...\n" 
        # Remove detected vector regions
        $DeConProg --fasta $DIR_sif/$rnatri.Sifted.fasta --csv $DIR_dec/$rnatri.decont --out $DIR_sif --keep && 
        mv -f $DIR_sif/$rnatri.Sifted.Sifted.fasta $DIR_sif/$rnatri.Sifted.fasta &&         
        echo -e "\nFinished. Clean assembly saved to $DIR_sif/$rnatri.Sifted.fasta...\n$(date)\n" 

        echo -e "\n========================\nScreening $DIR_sif/$rnatri.Sifted.fasta to confirm no vector sequences...\n$(date)\n========================\n" 
        # Re-screen to verify vector removal
        rm $DIR_dec/$rnatri.decont && 
        blastn -db $UniVecDB -query $DIR_sif/$rnatri.Sifted.fasta -out $DIR_dec/$rnatri.decont -outfmt "6 qseqid sseqid qstart qend staxids" -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 && 
        echo -e "\nFinished screening $DIR_sif/$rnatri.Sifted.fasta.\n$(date)\n" 
        
        ((counter++))
    done
    
    echo -e "\n########################\nCleaning Vectors from $rnatri is accomplished after $counter rounds.\n$(date)\n########################\n" 
}

# runFindAlien: Identify sequences from non-target organisms using DIAMOND blastx
# Searches clean assembly against NR protein database to find contaminating taxa.
# Input: $1 = full path to vector-cleaned FASTA file
runFindAlien(){

    local rnatri=$(echo $1 | sed -e "s#$DIR_sif\/##g;s#.Sifted.fasta##g")

    echo -e "\n########################\nSearching contaminated sequences...\n$(date)\n########################\n" 

    # Check if DIAMOND search already completed; if not, run it
    if [[ ! -f $DIR_dec/$rnatri.decontX ]] || [[ $(ls -l $DIR_dec/$rnatri.decontX | awk '{ print $5 }') -eq 0 ]]; then 
        echo -e "\n========================\nSearching unwanted transcripts by DIAMOND blastx against NR database....\n$(date)\n========================\n" && 
        diamond blastx --threads $DNP --db $NRDiaDB --query $1 --out $DIR_dec/$rnatri.decontX --outfmt 6 qseqid sseqid qstart qend score evalue pident ppos staxids --max-target-seqs 10 --fast && 
        echo -e "\nFinished screening $DIR_sif/$rnatri.Sifted.fasta.\n$(date)\n" 
    else
        echo -e "\n########################\nScreening $rnatri.Sifted.fasta has already been accomplished, no need to do it again.\n$(date)\n########################\n" 
    fi
}

# runDeCont: Remove sequences from non-target taxa using taxonomic filtering
# Uses NCBI taxonomy to filter out alien sequences, keeping only target organism taxa.
# Input: $1 = full path to vector-cleaned FASTA file
runDeCont(){

    local rnatri=$(echo $1 | sed -e "s#$DIR_sif\/##g;s#.Sifted.fasta##g")

    # Select decontamination program and validation rank based on organism type
    if [[ $EcdysFasta =~ $rnatri ]]; then 
        local DeConProg=/To/Your/Directory/nematoda/data/DeContaminRNA-tx.py
        local DeConVad=Ecdysozoa                        # Validate against Ecdysozoa for these organisms
    else
        local DeConProg=/To/Your/Directory/nematoda/data/DeContaminRNA-rk.py
        local DeConVad=phylum                           # Validate against organism's phylum for others
    fi
    
    echo -e "\n########################\nRemoval of contaminated sequences...\n$(date)\n########################\n" 

    echo -e "\n========================\nRemoving alien transcripts from the assembly $rnatri ...\n$(date)\n========================\n" && 
    $DeConProg --fasta $DIR_sif/$rnatri.Sifted.fasta --csv $DIR_dec/$rnatri.decontX --out $DIR_sif --sensitive 0 --valid $DeConVad --db $TaxDumpDB &&
    echo -e "\nFinished. Clean assembly saved to $DIR_sif/$rnatri.Sifted.DeContX.fasta...\n$(date)\n" 

    echo -e "\n########################\nCleaning alien transcripts from $rnatri is accomplished...\n$(date)\n########################\n" 
}

chmod -R 777 /To/Your/Directory/nematoda/data/*.py 

source $DIR_cnd/bin/activate $DIR_cnd/envs/blastenv 

# ============================================================================
# Stage 1: Vector Removal
# ============================================================================
# Setup file descriptor for process pool management
# Multithread preparation
tmp_fifofile="/tmp/$$.fifo" 
trap "exec 9>&-;exec 9<&-;exit 0" 2 
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$FNP;i++))
do
    echo >&9
done

# Vector removal processes
for rnatri in $InputFasta
do
    read -u9
    {

    # Skip if vector removal already accomplished
    if [[ ! -f $DIR_dec/$rnatri.decont ]] || [[ $(ls -l $DIR_dec/$rnatri.decont | awk '{ print $5 }') -gt 0 ]]; then 
        runDeVect $DIR_asm/$rnatri.fasta >> $DIR_dec/$rnatri.DecontamRNA.6.log 2>&1 
    else
        echo -e "\n########################\nCleaning Vectors from $rnatri has already been accomplished, no need to do it again.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log 
    fi 

    echo >&9 
    } & 
done
wait

# ============================================================================
# Stage 2: Alien Sequence Detection
# ============================================================================
# Multithread preparation
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$NT;i++))
do
    echo >&9
done

# Alien sequence detection processes
for rnatri in $InputFasta
do
    read -u9
    {

    # Detect non-target taxa using DIAMOND
    while [[ ! -f $DIR_dec/$rnatri.decontX ]] || [[ $(ls -l $DIR_dec/$rnatri.decontX | awk '{ print $5 }') -eq 0 ]] 
    do
        runFindAlien $DIR_sif/$rnatri.Sifted.fasta >> $DIR_dec/$rnatri.DecontamRNA.6.log 2>&1 
    done

    echo >&9 
    } & 
done
wait

source $DIR_cnd/bin/deactivate 

source $DIR_cnd/bin/activate $DIR_cnd/envs/rnaseqenv

# ============================================================================
# Stage 3: Taxonomic Filtering, ORF Prediction, and Visualization
# ============================================================================
# Multithread preparation
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$FNP;i++))
do
    echo >&9
done

# 
SiftFasta=$(ls $DIR_dec/*_trinity.Trinity.decontX | sed -e "s#.decontX##g;s#$DIR_dec\/##g")

# Decontamination, ORF prediction, and visualization processes
for rnatri in $SiftFasta
do
    read -u9
    {

    # Apply taxonomic filtering to remove non-target organisms
    while [[ ! -f $DIR_sif/$rnatri.Sifted.DeContX.fasta ]] || [[ $(ls -l $DIR_sif/$rnatri.Sifted.DeContX.fasta | awk '{ print $5 }') -lt 100 ]] 
    do
        runDeCont $DIR_sif/$rnatri.Sifted.fasta >> $DIR_dec/$rnatri.DecontamRNA.6.log 2>&1 
    done
    
    # Predict open reading frames (ORFs) and translate to proteins
    if [[ ! -f $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa ]]; then
        echo "\n========================\nSelecting the longest isoform.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log && 
        perl $DIR_cnd/envs/rnaseqenv/opt/trinity-2.13.2/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $DIR_sif/$rnatri.Sifted.DeContX.fasta > $DIR_sif/$rnatri.Sifted.DeContX.L.fasta && 
        
        echo "\n========================\nTransDecoding $rnatri (longest ORF search) starts...\n$(date)\n========================\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log && 
        TransDecoder.LongOrfs -t $DIR_sif/$rnatri.Sifted.DeContX.L.fasta >> $DIR_dec/$rnatri.DecontamRNA.6.log 2>&1 && 
        TransDecoder.Predict -t $DIR_sif/$rnatri.Sifted.DeContX.L.fasta >> $DIR_dec/$rnatri.DecontamRNA.6.log 2>&1 && 
        mv $DIR_cur/"$rnatri"*.transdecoder.* $DIR_tdc && 
        mv $DIR_tdc/"$rnatri"*.transdecoder.pep $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa && 
        rm -rf $DIR_cur/$rnatri*.transdecoder_dir* && 
        echo "\n========================\nTransDecoding $rnatri (longest ORF prediction) finished.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log 
    else
        echo -e "\n########################\nTransDecoding has already been accomplished, no need to do it again.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log 
    fi

    # Generate pie chart visualization of transcript taxonomic composition
    if [[ ! -f $DIR_Trc/$rnatri.png ]]; then
        echo -e "\n========================\nGenerating taxonomic composition visualization for $rnatri...\n$(date)\n========================\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log && 
        $SeqOrinPiePlot --infile $DIR_dec/$rnatri.decontX --out $DIR_Trc --db $TaxDumpDB >> $DIR_dec/$rnatri.DecontamRNA.6.log 2>&1 && 
        echo -e "\n========================\nVisualization complete for $rnatri.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log         
    else
        echo -e "\n########################\nVisualization has already been generated, no need to do it again.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log 
    fi

    echo -e "\n########################\nAll data processing for $rnatri has been completed successfully.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.DecontamRNA.6.log 

    echo >&9 
    } & 
done
wait


