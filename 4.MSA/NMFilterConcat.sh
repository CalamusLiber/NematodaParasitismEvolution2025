#!/bin/bash
#SBATCH --job-name=NMBU.x
#SBATCH --exclusive
#SBATCH --partition=mwvdk
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=350GB

# =============================================================================
# NMFilterConcat.sh - BUSCO Gene Tree Inference and Quality Filtering
# =============================================================================
#
# DESCRIPTION:
#   This SLURM script performs comprehensive BUSCO gene tree inference and quality
#   filtering for phylogenomic analysis. It processes both normally trimmed and
#   stationary-based trimmed alignments, performs likelihood mapping, symmetry
#   testing, and gene tree quality assessment using IQ-TREE and PhyKit.
#
# PIPELINE OVERVIEW:
#   1. IQ-TREE gene tree inference for normal trimming (nom)
#   2. IQ-TREE gene tree inference for stationary-based trimming (sbt)
#   3. Gene tree quality assessment using PhyKit metrics
#   4. Locus selection based on multiple quality criteria
#   5. Supermatrix generation for downstream phylogenetics
#
# INPUTS:
#   - $DIR_matrix/loci/nom/ : Normal trimmed BUSCO alignments
#   - $DIR_matrix/loci/sbt/ : Stationary-based trimmed alignments
#   - $DIR_matrix/LociFilterMeasure.*.list : Alignment quality metrics
#
# OUTPUTS:
#   - $DIR_matrix/genetree/nom/ : Gene trees from normal trimming
#   - $DIR_matrix/genetree/sbt/ : Gene trees from stationary-based trimming
#   - $DIR_matrix/GeneTreeMeasure.*.list : Gene tree quality metrics
#   - $DIR_matrix/AllBackbone*.rmsp : Filtered locus lists for supermatrix
#
# DEPENDENCIES:
#   - IQ-TREE2 (for tree inference and model selection)
#   - PhyKit (for tree quality metrics)
#   - Python scripts: LociSelection_aliba.py, LociSelection_treba.py
#
# AUTHOR: [Original author]
# MODIFIED: [Date] - Added comprehensive documentation
# =============================================================================

# Set bash options for robust error handling
set -euo pipefail

# Calculate parallel processing parameters from SLURM environment
NP=$SLURM_CPUS_PER_TASK  # Number of CPUs per task (8)
NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE]  # Total number of tasks (16)

# Load required software modules
module load apps/iqtree/2.1.3
# module load apps/bmge/1.12  # Commented out - BMGE loaded elsewhere if needed
# module load apps/bmge/1.12  # Duplicate commented module

# Define directory paths for data and tools
DIR_cnd=/To/Your/Directory/anaconda3          # Conda environment directory
DIR_united=/To/Your/Directory/nematoda/data/busco/united  # BUSCO results directory
DIR_msa=/To/Your/Directory/nematoda/msa       # MSA working directory
DIR_align=$DIR_msa/align                      # Alignment output directory
DIR_trim=$DIR_msa/trim                        # Trimmed alignment directory
DIR_matrix=$DIR_msa/matrix                    # Matrix and tree output directory
DIR_log=$DIR_msa/log                          # Log file directory

# Define tool paths
IQTREE=/sw/apps/iqtree-2.1.3-Linux/bin/iqtree2  # IQ-TREE2 executable path
LOCISELECTION_AB=$DIR_msa/LociSelection_aliba.py  # Alignment-based locus selection script
LOCISELECTION_TB=$DIR_msa/LociSelection_treba.py  # Tree-based locus selection script

# Create necessary output directories
mkdir -p $DIR_matrix/loci
mkdir -p $DIR_matrix/loci/nom  # Directory for normal trimmed loci
mkdir -p $DIR_matrix/loci/sbt  # Directory for stationary-based trimmed loci
mkdir -p $DIR_matrix/genetree
mkdir -p $DIR_matrix/genetree/nom
mkdir -p $DIR_matrix/genetree/sbt
mkdir -p $DIR_matrix/lmap

# Set permissions and activate conda environment
chmod -R 777 $LOCISELECTION_AB $LOCISELECTION_TB
source $DIR_cnd/bin/activate $DIR_cnd/envs/msaenv

# =============================================================================
# SECTION 1: Normal BMGE Trimmed Sequences Processing
# =============================================================================
# This section processes normally trimmed BUSCO alignments through:
# 1. Locus selection based on alignment quality metrics
# 2. Likelihood mapping to assess phylogenetic signal
# 3. Symmetry testing for model selection
# 4. Gene tree inference with bootstrap support

### Normal BMGE trimmed sequences
## Select Loci with PHYKIT results
echo -e "=======================\nSelect Loci for IQTREE-based alignment filtering...\nNormal trimming\n$(date)\n=======================\n\n" >> $DIR_log/MSALociSelection.log
$LOCISELECTION_AB -i $DIR_matrix/LociFilterMeasure.nom.list -o $DIR_matrix/RetainedLoci.nom.list

## IQTREE filtering
cd $DIR_matrix
RetainedLoci_nom=$(cat $DIR_matrix/RetainedLoci.nom.list | sed -e "s#$DIR_trim\/nom\/##g")

# Set up FIFO for parallel processing coordination
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile
rm $tmp_fifofile
for ((i=1;i<=$[$NT*$NP];i++))
do
    echo >&9
done

# Perform likelihood mapping if not already done
if [[ ! -f $DIR_matrix/Lmap50Loci.nom.list ]]; then
    for gene in $RetainedLoci_nom
    do
        read -u9  # Acquire semaphore token for parallel processing
        {
            # Log the start of likelihood mapping for this gene
            echo -e "=======================\nLikelihood mapping with IQTREE2 start...\n$gene - Normal trimming\n$(date)\n=======================\n\n" >> $DIR_log/MSALociSelection.log &&

            # Run IQ-TREE likelihood mapping if output doesn't exist
            # -s: input alignment file
            # -msub nuclear: test nuclear models only
            # -lmap 100000: perform likelihood mapping with 100k quartets
            # -n 0: no tree search, only likelihood mapping
            # --prefix: output file prefix
            [[ -f $DIR_matrix/lmap/$gene.iqtree ]] ||
            $IQTREE -s $DIR_trim/nom/$gene -msub nuclear -lmap 100000 -n 0 --prefix $DIR_matrix/lmap/$gene

            # Extract percentage of fully resolved quartets and filter loci with >=50%
            Pquartet=$(grep 'Number of fully resolved quartets (regions 1+2+3):' $DIR_matrix/lmap/$gene.iqtree | awk -F '[=%]' '{print $2}') &&
            [[ $(echo "$Pquartet >= 50" | bc) -eq 1 ]] && echo -e "$DIR_trim/nom/$gene\t$Pquartet" >> $DIR_matrix/Lmap50Loci.nom.list &&

            # Log completion of likelihood mapping
            echo -e "=======================\n$gene - Normal trimming\nLikelihood mapping with IQTREE2 accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/MSALociSelection.log 
            
        echo >&9
        } &
    done
    wait
fi

## IQTREE symmetry testing and filtering
cd $DIR_matrix &&

# Perform model selection and symmetry testing if not already done
if [[ ! -f $DIR_matrix/loci/nom.best_scheme.nex ]]; then
    # Copy high-quality loci (from likelihood mapping) to loci directory
    cp $(cat $DIR_matrix/Lmap50Loci.nom.list | awk '{print $1}') $DIR_matrix/loci/nom &&
    # Run IQ-TREE partition analysis with symmetry testing
    # -p: partition file mode
    # -msub nuclear: test nuclear models only
    # --symtest-remove-bad: perform symmetry test and remove poorly fitting partitions
    # -T: number of threads
    $IQTREE -p $DIR_matrix/loci/nom -msub nuclear --symtest-remove-bad -T $[$NT*$NP/2]
fi

## IQTREE tree-building for each locus (gene trees)
# Extract list of loci that passed symmetry testing
IQT_list=$(grep '  charset ' $DIR_matrix/loci/nom.best_scheme.nex | awk '{print $2}' | sed -e "s#.mft.bmge.nom.faa##g")

# Set up FIFO for parallel gene tree inference
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile
rm $tmp_fifofile
for ((i=1;i<=$NT;i++))
do
    echo >&9
done

# Build gene trees for each locus that passed filtering
for gene in $IQT_list
do
    read -u9  # Acquire semaphore token
    {
        if [[ ! -f $DIR_matrix/genetree/nom/$gene.contree ]]; then
            # Extract the best-fit model for this locus from the partition scheme
            submdl=$(grep ": $gene\." $DIR_matrix/loci/nom.best_scheme.nex | awk '{print $1}' | sed 's/://g')
            # Run IQ-TREE tree inference with bootstrap support
            # -s: input alignment
            # -m: substitution model
            # -B 1000: 1000 ultrafast bootstrap replicates
            # -T: number of threads
            # --prefix: output file prefix
            $IQTREE -s $DIR_matrix/loci/nom/$gene.mft.bmge.nom.faa --prefix $DIR_matrix/genetree/nom/$gene -m $submdl -B 1000 -T $NP
        fi
        echo >&9  # Release semaphore token
    } &
done
wait

# =============================================================================
# SECTION 2: Stationary-based BMGE Trimmed Sequences Processing
# =============================================================================
# This section processes stationary-based trimmed BUSCO alignments through:
# 1. Locus selection based on alignment quality metrics
# 2. Likelihood mapping to assess phylogenetic signal
# 3. Symmetry testing for model selection
# 4. Gene tree inference with bootstrap support

### Stationary-based BMGE trimmed sequences
## Select Loci with PHYKIT results
echo -e "=======================\nSelect Loci for IQTREE-based alignment filtering...\nStationary-based trimming\n$(date)\n=======================\n\n" >> $DIR_log/MSALociSelection.log
$LOCISELECTION_AB -i $DIR_matrix/LociFilterMeasure.sbt.list -o $DIR_matrix/RetainedLoci.sbt.list

## IQTREE filtering
cd $DIR_matrix
RetainedLoci_sbt=$(cat $DIR_matrix/RetainedLoci.sbt.list | sed -e "s#$DIR_trim\/sbt\/##g")

# Set up FIFO for parallel likelihood mapping
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile
rm $tmp_fifofile
for ((i=1;i<=$[$NT*$NP];i++))
do
    echo >&9
done

# Perform likelihood mapping if not already done
if [[ ! -f $DIR_matrix/Lmap50Loci.sbt.list ]]; then
    for gene in $RetainedLoci_sbt
    do
        read -u9  # Acquire semaphore token for parallel processing
        {
            # Log the start of likelihood mapping for this gene
            echo -e "=======================\nLikelihood mapping with IQTREE2 start...\n$gene - Stationary-based trimming\n$(date)\n=======================\n\n" >> $DIR_log/MSALociSelection.log &&

            # Run IQ-TREE likelihood mapping if output doesn't exist
            # -s: input alignment file
            # -msub nuclear: test nuclear models only
            # -lmap 100000: perform likelihood mapping with 100k quartets
            # -n 0: no tree search, only likelihood mapping
            # --prefix: output file prefix
            [[ -f $DIR_matrix/lmap/$gene.iqtree ]] ||
            $IQTREE -s $DIR_trim/sbt/$gene -msub nuclear -lmap 100000 -n 0 --prefix $DIR_matrix/lmap/$gene

            # Extract percentage of fully resolved quartets and filter loci with >=50%
            Pquartet=$(grep 'Number of fully resolved  quartets (regions 1+2+3):' $DIR_matrix/lmap/$gene.iqtree | awk -F '[=%]' '{print $2}') &&
            [[ $(echo "$Pquartet >= 50" | bc) -eq 1 ]] && echo -e "$DIR_trim/sbt/$gene\t$Pquartet" >> $DIR_matrix/Lmap50Loci.sbt.list &&

            # Log completion of likelihood mapping
            echo -e "=======================\n$gene - Stationary-based trimming\nLikelihood mapping with IQTREE2 accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/MSALociSelection.log

            echo >&9  # Release semaphore token
        } &
    done
    wait
fi

## IQTREE symmetry testing and filtering
cd $DIR_matrix && 

# Perform model selection and symmetry testing if not already done
if [[ ! -f $DIR_matrix/loci/sbt.best_scheme.nex ]]; then
    # Copy high-quality loci (from likelihood mapping) to loci directory
    cp $(cat $DIR_matrix/Lmap50Loci.sbt.list | awk '{print $1}') $DIR_matrix/loci/sbt &&
    # Run IQ-TREE partition analysis with symmetry testing
    # -p: partition file mode
    # -msub nuclear: test nuclear models only
    # --symtest-remove-bad: perform symmetry test and remove poorly fitting partitions
    # -T: number of threads
    $IQTREE -p $DIR_matrix/loci/sbt -msub nuclear --symtest-remove-bad -T $[$NT*$NP/2]
fi

## IQTREE tree-building for each locus (gene trees)
# Extract list of loci that passed symmetry testing
IQT_list=$(grep '  charset ' $DIR_matrix/loci/sbt.best_scheme.nex | awk '{print $2}' | sed -e "s#.mft.bmge.sbt.faa##g")

# Set up FIFO for parallel gene tree inference
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile
rm $tmp_fifofile
for ((i=1;i<=$NT;i++))
do
    echo >&9
done

# Build gene trees for each locus that passed filtering
for gene in $IQT_list
do
    read -u9  # Acquire semaphore token
    {
        if [[ ! -f $DIR_matrix/genetree/sbt/$gene.contree ]]; then
            # Extract the best-fit model for this locus from the partition scheme
            submdl=$(grep ": $gene\." $DIR_matrix/loci/sbt.best_scheme.nex | awk '{print $1}' | sed 's/://g')
            # Run IQ-TREE tree inference with bootstrap support
            # -s: input alignment
            # -m: substitution model
            # -B 1000: 1000 ultrafast bootstrap replicates
            # -T: number of threads
            # --prefix: output file prefix
            $IQTREE -s $DIR_matrix/loci/sbt/$gene.mft.bmge.sbt.faa --prefix $DIR_matrix/genetree/sbt/$gene -m $submdl -B 1000 -T $NP
        fi
        echo >&9  # Release semaphore token
    } &
done
wait

# =============================================================================
# SECTION 3: Gene Tree Quality Assessment and Supermatrix Generation
# =============================================================================
# This section assesses gene tree quality using PhyKit metrics and generates
# final supermatrices for phylogenomic analysis by filtering out problematic loci.

### Tests using gene trees
cd $DIR_matrix
echo "" > $DIR_matrix/nullroot  # Create empty root file for molecular clock analysis

gene_tree_nom=$(ls $DIR_matrix/genetree/nom/*.treefile)
gene_tree_sbt=$(ls $DIR_matrix/genetree/sbt/*.treefile)

# Set up FIFO for parallel gene tree quality assessment
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile
rm $tmp_fifofile
for ((i=1;i<=$[$NT*2];i++))
do
    echo >&9
done

# Assess quality metrics for normal trimming gene trees
for genetree in $gene_tree_nom
do
    read -u9  # Acquire semaphore token
    {
        # Calculate various tree quality metrics using PhyKit
        avgboots=$(echo $(phykit bipartition_support_stats $genetree) | awk '{print $2}') &&  # Average bootstrap support
        treeness=$(phykit treeness $genetree) &&  # Tree balance metric
        rcv=$(phykit relative_composition_variability $(echo $genetree | sed -e "s#$DIR_matrix\/genetree\/nom\/#$DIR_matrix\/loci\/nom\/#g;s#.treefile#.mft.bmge.nom.faa#g")) &&  # Compositional heterogeneity
        treenessdercv=$(echo "scale=4; $treeness / $rcv" | bc) &&  # Treeness/RCV ratio
        dvmc=$(phykit degree_of_violation_of_a_molecular_clock -t $genetree -r $DIR_matrix/nullroot) &&  # Molecular clock violation
        echo -e "$genetree,$avgboots,$treeness,$rcv,$treenessdercv,$dvmc" >> $DIR_matrix/GeneTreeMeasure.nom.list

        # Detect spurious sequences (potential paralogs or contaminants)
        sprseq=$(phykit spurious_seq $genetree -f 20 | sed "s/\t/,/g") &&
        if [[ $sprseq != None ]]; then
            for sp in $sprseq
            do
                echo -e $(echo $genetree | sed -e "s#$DIR_matrix\/genetree\/nom\/#$DIR_matrix\/loci\/nom\/#g;s#.treefile#.mft.bmge.nom.faa#g"),$sp >> $DIR_matrix/SpuriousHomolog.nom.list
            done
        fi

        echo >&9  # Release semaphore token
    } &
done
wait 

wait

# Assess quality metrics for stationary-based trimming gene trees
for genetree in $gene_tree_sbt
do
    read -u9  # Acquire semaphore token
    {
        # Calculate various tree quality metrics using PhyKit
        avgboots=$(echo $(phykit bipartition_support_stats $genetree) | awk '{print $2}') &&  # Average bootstrap support
        treeness=$(phykit treeness $genetree) &&  # Tree balance metric
        rcv=$(phykit relative_composition_variability $(echo $genetree | sed -e "s#$DIR_matrix\/genetree\/sbt\/#$DIR_matrix\/loci\/sbt\/#g;s#.treefile#.mft.bmge.sbt.faa#g")) &&  # Compositional heterogeneity
        treenessdercv=$(echo "scale=4; $treeness / $rcv" | bc) &&  # Treeness/RCV ratio
        dvmc=$(phykit degree_of_violation_of_a_molecular_clock -t $genetree -r $DIR_matrix/nullroot) &&  # Molecular clock violation
        echo -e "$genetree,$avgboots,$treeness,$rcv,$treenessdercv,$dvmc" >> $DIR_matrix/GeneTreeMeasure.sbt.list

        # Detect spurious sequences (potential paralogs or contaminants)
        sprseq=$(phykit spurious_seq $genetree -f 20 | sed "s/\t/,/g") &&
        if [[ $sprseq != None ]]; then
            for sp in $sprseq
            do
                echo -e $(echo $genetree | sed -e "s#$DIR_matrix\/genetree\/sbt\/#$DIR_matrix\/loci\/sbt\/#g;s#.treefile#.mft.bmge.sbt.faa#g"),$sp >> $DIR_matrix/SpuriousHomolog.sbt.list
            done
        fi

        echo >&9  # Release semaphore token
    } &
done
wait

## Generating new data set super matrix
cd $DIR_matrix

# Generate final filtered supermatrices using tree-based locus selection
# This integrates alignment quality metrics, gene tree quality metrics, and spurious sequence detection

# Generate supermatrix for normal trimming (full dataset with outgroups)
$LOCISELECTION_TB -a $DIR_matrix/LociFilterMeasure.nom.list -g $DIR_matrix/GeneTreeMeasure.nom.list -s $DIR_matrix/SpuriousHomolog.nom.list -o $DIR_matrix/TreeBasedDropLoci.nom.rmsp &&

# Generate supermatrix for Windows analysis (58 taxa including outgroups)
# -a: alignment quality metrics file
# -g: gene tree quality metrics file
# -s: spurious homologs file
# -o: output file for retained loci
# -e: taxa to exclude from supermatrix (outgroups for backbone analysis)
python3 $LOCISELECTION_TB -a $DIR_matrix/LociFilterMeasure.nom.copy.list -g $DIR_matrix/GeneTreeMeasure.nom.copy.list -s $DIR_matrix/SpuriousHomolog.nom.list -o $DIR_matrix/AllBackbone.nom.rmsp -e Acanthocheilonema_spirocauda,Acanthocheilonema_viteae,Allodiplogaster_sudhausi,Allodorylaimus_sp,Anaplectus_granulosus,Ancylostoma_caninum,Ancylostoma_ceylanicum,Ancylostoma_duodenale,Angiostrongylus_costaricensis,Anisakis_simplex,Aphelenchoides_besseyi,Aphelenchoides_fragariae,Ascaridia_galli,Ascaris_sp,Ascaris_suum,Bathylaimus_sp,Baylisascaris_schroederi,Brugia_malayi,Brugia_timori,Bursaphelenchus_mucronatus,Bursaphelenchus_xylophilus,Caenorhabditis_briggsae,Caenorhabditis_elegans,Ceramonema_sp,Cooperia_punctata,Crassicauda_magna,Cylicolaimus_sp,Cylicostephanus_goldi,Deladenus_siricidicola,Dictyocaulus_viviparus,Dipetalonema_caudispina,Diplogasteroides_magnus,Diploscapter_coronatus,Diploscapter_pachys,Dirofilaria_immitis,Ditylenchus_destructor,Ditylenchus_dipsaci,Dorylaimopsis_sp,Elaeophora_elaphi,Euchromadora_sp,Gammarinema_sp,Globodera_ellingtonae,Globodera_pallida,Gongylonema_pulchrum,Haemonchus_placei,Halicephalobus_mephisto,Halicephalobus_sp,Halichoanolaimus_dolichurus,Halomonhystera_disjuncta,Heligmosomoides_polygyrus,Heterodera_glycines,Heterodera_schachtii,Heterorhabditis_bacteriophora,Hoplolaimus_columbus,Koerneria_luziae,Laxus_oneistus,Leptonemella_sp,Levipalatum_texanum,Litomosoides_brasiliensis,Litomosoides_sigmodontis,Longidorus_elongatus,Madathamugadia_hiepei,Meloidogyne_enterolobii,Meloidogyne_hapla,Micoletzkya_masseyi,Monhystera_sp,Necator_americanus,Neocamacolaimus_parasiticus,Nippostrongylus_brasiliensis,Odontophora_sp,Oesophagostomum_dentatum,Onchocerca_flexuosa,Onchocerca_ochengi,Oncholaimidae_sp,Oscheius_sp,Panagrellus_redivivus,Panagrolaimus_davidi,Panagrolaimus_sp,Panagrolaimus_superbus,Paradraconema_sp,Parapristionchus_giblindavisi,Parascaris_equorum,Parascaris_univalens,Parastrongyloides_trichosuri,Parelaphostrongylus_tenuis,Pomponema_sp,Prionchulus_punctatus,Pristionchus_exspectatus,Pristionchus_mayeri,Pristionchus_pacificus,Propanagrolaimus_sp,Radopholus_similis,Rhabditophanes_sp,Rhigonematomorpha_sp,Rotylenchulus_reniformis,Setaria_digitata,Setaria_equina,Siphonolaimus_sp,Soboliphyme_baturini,Sphaerolaimus_sp,Steinernema_carpocapsae,Steinernema_scapterisci,Stephanolaimus_elegans,Stilbonema_sp,Strongyloides_papillosus,Strongyloides_stercoralis,Strongylus_vulgaris,Subanguina_moxae,Symplocostoma_sp,Syphacia_muris,Teladorsagia_circumcincta,Thelastomatoidea_sp,Thelazia_callipaeda,Thoonchus_sp,Thoracostomopsidae_sp,Toxascaris_leonina,Toxocara_canis,Trichinella_nativa,Trichinella_pseudospiralis,Trichostrongylus_colubriformis,Trichuris_muris,Trichuris_suis,Trichuris_trichiura,Turbatrix_aceti,Wuchereria_bancrofti,Xiphinema_index 

# run on BluePebble (43 taxa, only Nematoda)
$LOCISELECTION_TB -a $DIR_matrix/LociFilterMeasure.nom.list -g $DIR_matrix/GeneTreeMeasure.nom.list -s $DIR_matrix/SpuriousHomolog.nom.list -o $DIR_matrix/AllBackbone4.nom.rmsp -e PF_Amphimedon_queenslandica,CN_Nematostella_vectensis,ED_Anneissia_japonica,AN_Helobdella_robusta,PH_Echinococcus_granulosus,PH_Oisthorchis_felineus,ML_Acanthochitona_rubrolineata,ML_Mercenaria_mercenaria,ML_Octopus_sinensis,KR_Pycnophyes_sp,PA_Meiopriapulus_fijiensis,PA_Priapulus_caudatus,OP_Epiperipatus_trinidadensis,AP_Meridionale_flava,AP_Stegodyphus_dumicola,AP_Eriocheir_sinensis,AP_Daphnia_magna,AP_Megajapyx_major,AP_Sinella_curviseta,AP_Machilis_hrabei,AP_Timema_cristinae,AP_Tribolium_castaneum,NM_Nectonema_munidae,NM_Gordius_sp,NM_Paragordius_varius,TG_Richtersius_coronifer,TG_Hypsibius_dujardini,TG_Ramazzottius_varieornatus,Acanthocheilonema_spirocauda,Acanthocheilonema_viteae,Allodiplogaster_sudhausi,Allodorylaimus_sp,Anaplectus_granulosus,Ancylostoma_caninum,Ancylostoma_duodenale,Angiostrongylus_costaricensis,Anisakis_simplex,Aphelenchoides_besseyi,Aphelenchoides_fragariae,Ascaris_sp,Ascaris_suum,Bathylaimus_sp,Baylisascaris_schroederi,Brugia_malayi,Brugia_timori,Bursaphelenchus_mucronatus,Bursaphelenchus_xylophilus,Caenorhabditis_briggsae,Caenorhabditis_elegans,Cooperia_punctata,Cylicostephanus_goldi,Deladenus_siricidicola,Dictyocaulus_viviparus,Dipetalonema_caudispina,Diplogasteroides_magnus,Diploscapter_coronatus,Dirofilaria_immitis,Ditylenchus_destructor,Ditylenchus_dipsaci,Elaeophora_elaphi,Euchromadora_sp,Gammarinema_sp,Globodera_ellingtonae,Globodera_pallida,Haemonchus_placei,Halicephalobus_mephisto,Halicephalobus_sp,Halichoanolaimus_dolichurus,Halomonhystera_disjuncta,Heligmosomoides_polygyrus,Heterodera_glycines,Heterodera_schachtii,Hoplolaimus_columbus,Koerneria_luziae,Laxus_oneistus,Leptonemella_sp,Levipalatum_texanum,Litomosoides_brasiliensis,Litomosoides_sigmodontis,Longidorus_elongatus,Madathamugadia_hiepei,Meloidogyne_enterolobii,Meloidogyne_hapla,Micoletzkya_masseyi,Monhystera_sp,Necator_americanus,Neocamacolaimus_parasiticus,Nippostrongylus_brasiliensis,Oesophagostomum_dentatum,Onchocerca_flexuosa,Onchocerca_ochengi,Oncholaimidae_sp,Oscheius_sp,Panagrellus_redivivus,Panagrolaimus_davidi,Panagrolaimus_sp,Panagrolaimus_superbus,Paradraconema_sp,Parapristionchus_giblindavisi,Parascaris_equorum,Parascaris_univalens,Parastrongyloides_trichosuri,Parelaphostrongylus_tenuis,Pomponema_sp,Prionchulus_punctatus,Pristionchus_exspectatus,Pristionchus_mayeri,Pristionchus_pacificus,Propanagrolaimus_sp,Radopholus_similis,Rhabditophanes_sp,Rotylenchulus_reniformis,Setaria_equina,Siphonolaimus_sp,Soboliphyme_baturini,Sphaerolaimus_sp,Steinernema_carpocapsae,Steinernema_scapterisci,Stephanolaimus_elegans,Stilbonema_sp,Strongyloides_papillosus,Strongyloides_stercoralis,Strongylus_vulgaris,Subanguina_moxae,Symplocostoma_sp,Syphacia_muris,Teladorsagia_circumcincta,Thelastomatoidea_sp,Thoonchus_sp,Thoracostomopsidae_sp,Toxascaris_leonina,Toxocara_canis,Trichinella_nativa,Trichinella_pseudospiralis,Trichostrongylus_colubriformis,Trichuris_muris,Trichuris_suis,Trichuris_trichiura,Turbatrix_aceti,Wuchereria_bancrofti,Xiphinema_index 

echo -e "=======================\nAll IQTREE-based alignment filtering jobs are accomplished.\n$(date)\n=======================\n\n" >> $DIR_log/MSALociSelection.log 
