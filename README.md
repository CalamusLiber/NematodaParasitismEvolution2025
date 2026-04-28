# Nematoda Phylogenomics & Evolutionary Analysis Pipeline (Draft-20260428)

## Overview

This repository accompanies our study on the deep evolutionary history of Nematoda (roundworms). It contains the full bioinformatics workflow we used to reconstruct nematode phylogeny, estimate divergence times, and trace the ancestral origins of parasitic life history — starting from raw sequencing reads and ending with publication-ready figures.

We built the pipeline to be modular: each stage can be run independently once its inputs are in place, and most computationally heavy steps are designed to run on an HPC cluster via SLURM or PBS. If you are adapting this for your own data, you will mostly need to update the file paths and sample names scattered through the scripts.

The workflow runs through seven sequential stages:

```
0.Environment  →  1.SequenceAcquisition  →  2.Assembly  →  3.Prediction
                                                                    ↓
7.AncestralStateReconstruction  ←  6.MolecularClock  ←  5.Phylogenetics  ←  4.MSA
```

When utilizing the methods, algorithms, or codes developed for this study, even if they are not related to nematodes, please be sure to cite our paper:

```
Lü, L., K. De Baets, M. Giacomelli, M. E. Rossi, O. Holovachov, D. Pisani, P. C. J. Donoghue (2026) The timescale of nematode evolution reveals late and convergent radiations of parasitism. (in press).
```

---

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Installation Guide](#installation-guide)
3. [Pipeline Overview & Script Descriptions](#pipeline-overview--script-descriptions)
4. [Demo](#demo)
5. [Instructions for Use](#instructions-for-use)
6. [Expected Outputs](#expected-outputs)
7. [Expected Run Times](#expected-run-times)
8. [Dependencies Summary](#dependencies-summary)

---

## System Requirements

### Hardware

Most steps in this pipeline genuinely need an HPC cluster — genome assembly, gene prediction, tree inference, and MCMCtree are simply not practical on a laptop. That said, the R and Python analysis scripts (Stages 6–7) run fine on a well-equipped desktop.

- **HPC cluster** with SLURM or PBS/Torque (required for Stages 2–5)
- Compute nodes: ≥ 24 CPU cores, ≥ 100 GB RAM per node; genome assembly jobs need ≥ 256 GB RAM
- Storage: ≥ 2 TB free space for the full dataset (raw reads, assemblies, alignments, and MCMC output accumulate fast)

### Operating System

- **Linux** — all scripts were developed and tested on CentOS 7 and Ubuntu 20.04+ on HPC; this is the expected environment.
- **Windows** — the R scripts (Stages 6–7) can be run interactively on a Windows desktop; Python scripts have not been tested on Windows and may need minor path adjustments.
- **Mac** — none of the code has been tested on macOS; it may work for the R/Python components but is untested.

### Core Software
- **Conda / Anaconda** ≥ 23.x — we use Conda to manage all tool environments, which avoids most dependency conflicts.
- **Python** ≥ 3.9 (3.10 or 3.12 recommended) for all custom Python scripts. Notably, certain third-party shoftware may depend on Python 2.7.
- **R** ≥ 4.2 for all the post-time-tree-estimation analyses.

---

## Guide for Environment Configuration

Our pipeline is not a standalone application or an integrated software package. The scripts provided here serve as examples to demonstrate how we implemented the methods used and developed during our study. Therefore, we do not provide a one-click installer or a standard configure/make installation process.

### 1. Set Up the Project Directory

Open `0.Environment/init.sh` and set `BaseinDir` to wherever you want to store the analysis data. Then run:

```bash
bash 0.Environment/init.sh
```

This creates the full directory tree (`data/`, `msa/`, `phylo/`, `timing/`) in one go, so you don't have to make folders manually as you go.

### 2. Install All Conda Environments

All external tools are managed through Conda. Running the installation script takes 30–60 minutes depending on your connection, but you only need to do it once:

```bash
bash 0.Environment/CondaInstallation.sh
```

This creates the following environments:

| Environment | Purpose |
|---|---|
| `blastenv` | BLAST+ and DIAMOND sequence search |
| `makerenv` | MAKER gene prediction pipeline |
| `repmaskenv` | RepeatMasker repeat annotation |
| `braker2env` | BRAKER2 automated gene prediction |
| `buscoenv` | BUSCO genome completeness |
| `orthofinderenv` | OrthoFinder orthology detection |
| `assemblingenv` | SPAdes / SOAPdenovo2 genome assembly |
| `rnaseqenv` | Trinity / SOAPdenovo-Trans transcriptome |
| `seqqc` | FastP / FastQC quality control |
| `msaenv` | MAFFT / BMGE alignment & trimming |
| `newpy` | Python 3.12 general purpose |

### 3. Install Core R Packages

Most R packages can be installed directly from CRAN in a single call. Open R (or RStudio) and run:

```r
install.packages(c(
  # Phylogenetics & tree handling
  "ape", "geiger", "phytools", "ips",
  # Visualisation
  "ggplot2", "patchwork", "ggalt", "grid",
  # Diversification & ancestral state
  "hisse",
  # Statistical modelling
  "vegan", "car", "agricolae", "plyr", "dplyr",
  "sn", "pracma",
  # Machine learning
  "xgboost", "rpart", "caret",
  # Tidymodels ecosystem
  "tidymodels", "recipes", "parsnip", "brulee",
  # Deep learning & clustering
  "torch", "dbscan", "Rtsne", "fpc",
  # Parallelism
  "doParallel", "future",
  # Utilities
  "readxl"
))
```

Two packages require **Bioconductor** (install `BiocManager` first if you don't have it):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
```

Two packages are best installed from **GitHub** for the most up-to-date versions:

```r
# devtools or remotes is required
install.packages("remotes")

# HiSSE — hidden-state speciation & extinction
remotes::install_github("thej022214/hisse")

# tfprobability — TensorFlow Probability interface for R
# Requires a working TensorFlow installation; see https://rstudio.github.io/tfprobability/
install.packages("tfprobability")
# or
remotes::install_github("rstudio/tfprobability")
# then install TensorFlow and TensorFlow Probability python modules
library(tfprobability)
install_tfprobability()
```

> **Note:** `torch` and `brulee` require a one-time backend download after installation:
> ```r
> library(torch); torch::install_torch()
> ```
> `tfprobability` requires Python and TensorFlow to be available (managed via `reticulate`).

---

### 4. Build BLAST / DIAMOND Databases

Before you can run decontamination or annotation, you need local copies of the NCBI `nr` protein database and `UniVec.fasta`. Once downloaded, build the indexed databases:

```bash
bash 0.Environment/BuildingBlast+databases.sh
```

This builds BLAST+ (`nr`, `UniVec`) and DIAMOND (`nr`, `uniref100`) databases. These are large — `nr` alone is several hundred GB — so make sure you have the storage before starting.

---

## Pipeline Overview & Script Descriptions

Below is a stage-by-stage breakdown of every script in the repository, what it does, and what it depends on. Shell scripts marked with a scheduler (SLURM / PBS) need to be submitted to a cluster; the rest can be run directly.

---

### Stage 0 — Environment Setup (`0.Environment/`)

Three simple setup scripts that prepare your workspace and install everything else.

| Script | Language | Function |
|---|---|---|
| `init.sh` | Shell | Creates the project directory hierarchy |
| `CondaInstallation.sh` | Shell | Installs all Conda environments and tools |
| `BuildingBlast+databases.sh` | Shell | Builds BLAST+ and DIAMOND sequence databases |

---

### Stage 1 — Sequence Acquisition (`1.SequenceAcquisition/`)

These scripts search NCBI for available nematode genome assemblies and RNA-seq data, and then download them in parallel. The two Python scripts are designed to be run interactively (or from the command line) to generate the download lists; the shell script then handles the actual downloading on the cluster.

| Script | Language | Function |
|---|---|---|
| `get_assembly_info_20210819.py` | Python | Queries NCBI Assembly via BioPython Entrez; retrieves accessions, metadata, and FTP download URLs; outputs to Excel |
| `get_sra_info_max20210831.py` | Python | Queries NCBI SRA; fetches paired-end RNA-seq and WGS run metadata; filters by library strategy; outputs to Excel |
| `GenomeDownload.sh` | Shell (SBATCH) | Downloads genome assemblies and SRA data in parallel (up to 20 threads); can resume from the last processed row if interrupted |

**Key Python dependencies:** `biopython`, `beautifulsoup4`, `pandas`, `tqdm`

---

### Stage 2 — Assembly (`2.Assembly/`)

This stage covers both genome (WGS) and transcriptome (RNA-seq) assembly, plus the decontamination steps that follow. Decontamination is important here because many nematode samples carry bacterial or fungal reads that can end up in assemblies if left unchecked.

#### Shell scripts

| Script | Scheduler | Function |
|---|---|---|
| `Nematoda_AfterQC.4.3.sh` | — | Batch RNA-seq quality control using AfterQC |
| `Nematoda_spades.0.pbs` | PBS | WGS genome assembly pipeline: deduplication → quality trimming → depth normalization (BBtools) → SPAdes assembly |
| `NMSOAP_trans.sh` / `NMSOAP_trans.x.sh` | — | Transcriptome assembly using SOAPdenovo-Trans |
| `NMTrinity.9.sh` / `NMTrinity.9.x.sh` | SLURM | Per-sample Trinity de novo transcriptome assembly, AfterQC trimming, and TransDecoder ORF prediction |
| `NMTrinityQuality.sh` / `NMTransAssemblyQuality.sh` | — | Post-assembly quality statistics |
| `Trinity.DeConX.longestORF.sh` | — | Extracts the longest ORF per transcript from decontaminated Trinity assemblies |
| `UniVec.Decontamination.6.sh` | — | Screens for vector contamination using BLAST against the NCBI UniVec database |
| `soap2trinity.sh` | — | Converts SOAPdenovo-Trans output FASTA headers to Trinity-compatible naming format |
| `NMDrawSeqOri.slm` | SLURM | Submits sequence-origin analysis jobs |

#### Python scripts

| Script | Function |
|---|---|
| `DeContaminRNA-rk.py` | Removes non-nematode contigs from RNA assemblies using rank-based taxonomic filtering against DIAMOND/BLAST hits and the NCBI taxonomy tree |
| `DeContaminRNA-tx.py` | Alternate decontamination implementation using taxon-name-based filtering (useful when rank resolution is ambiguous) |
| `RmVecRNA.mp.py` | Multi-threaded removal of BLAST-identified vector regions from assembled contigs, trimming only the contaminated segments and preserving the clean flanks |
| `SeqOrigin.py` | Analyses taxonomic composition of assembled transcripts and classifies sequences by phylum |
| `SeqOriginPiePlot.py` | Generates pie charts of taxonomic composition; supports both live NCBI Taxonomy queries and a local taxonomy database for offline use |

**Key Python dependencies:** `biopython`, `pandas`, `matplotlib`, `numpy`, `multiprocessing`

---

### Stage 3 — Gene Prediction (`3.Prediction/`)

Here we mask repetitive elements, assess genome completeness with BUSCO, and run gene prediction using MAKER (multi-round, integrating RNA-seq evidence with Augustus and SNAP). The `UniteBUSCOs.py` script then collects the per-species BUSCO results into the unified gene sets that feed into Stage 4.

| Script | Scheduler | Function |
|---|---|---|
| `NMRepeatMasker.sh` | — | Repeat annotation with RepeatMasker and RepeatModeler |
| `GffHardMasker.sh` | — | Shell-based hard masking of genome FASTA sequences using RepeatMasker GFF coordinates (replaces masked intervals with `N`). It runs at a slow crawl |
| `NMBusco.sh` | SLURM | BUSCO genome completeness assessment |
| `NMmetaeuk.sh` / `NMmetaeuk.1.sh` | — | Gene prediction using MetaEuk against the UniRef90 protein database |
| `NMPredMaker.6.a.sh` | SLURM | Dispatcher: submits one MAKER job per genome, allocating resources adaptively based on genome size (>200 MB → 32 cores / 350 GB RAM; smaller genomes → 24 cores / 100 GB RAM) |
| `NMPredMaker.6.a.x.sh` | SLURM | Worker: runs the full multi-round MAKER pipeline (BLAST evidence integration → Augustus → SNAP) for a single genome |

#### Python scripts

| Script | Function |
|---|---|
| `GffRepeatMasker.py` | Applies repeat masking to FASTA sequences based on RepeatMasker GFF coordinates; choose between hard masking (replace with N) or soft masking (lowercase) |
| `UniteBUSCOs.py` | Consolidates BUSCO protein predictions across all taxa — selects the longest isoform per locus per taxon, writes unified FASTA files, and produces a presence/absence heatmap across all species |

**Key Python dependencies:** `biopython`, `pandas`, `matplotlib`, `numpy`

---

### Stage 4 — Multiple Sequence Alignments (`4.MSA/`)

Every BUSCO locus gets aligned independently with MAFFT, trimmed with BMGE (we run two trimming strategies in parallel to assess sensitivity), and then evaluated with PhyKit. The two Python loci-selection scripts implement different filtering philosophies — use `LociSelection_aliba.py` for a quick, metric-based filter, or `LociSelection_treba.py` for a more thorough approach that also uses gene-tree evidence to catch spurious homologs.

| Script | Scheduler | Function |
|---|---|---|
| `NMAlignTrim.sh` | SLURM (4 nodes × 4 tasks × 8 CPUs) | Parallel MAFFT alignment (E-INS-i) + BMGE trimming (two strategies) + PhyKit quality metrics for all BUSCO loci |
| `NMFilterConcat.sh` | — | Filters aligned loci and concatenates the selected set into a supermatrix |

#### Python scripts

| Script | Function |
|---|---|
| `LociSelection_aliba.py` | Quick loci filtering using quantile thresholds on alignment quality metrics (taxon count, parsimony-informative sites, alignment length, RCV, error rate) |
| `LociSelection_treba.py` | Advanced filtering combining alignment metrics with gene-tree quality scores (bootstrap, treeness, DVMC); removes spurious homologs identified in gene trees; outputs a concatenated PHYLIP/FASTA supermatrix with partition files |

**Key Python dependencies:** `biopython`, `pandas`, `numpy`, `matplotlib`

---

### Stage 5 — Phylogenetics (`5.Phylogenetics/`)

We inferred trees using two complementary methods: maximum likelihood with IQ-TREE 2, and Bayesian inference with PhyloBayes MPI (CAT-GTR model). The PhyloBayes runs are slow — expect days to weeks for large datasets — but the convergence diagnostics (`NMPBcomp.sh`) make it straightforward to check whether chains have mixed.

#### IQ-TREE (`IQTREE/`)

| Script | Scheduler | Function |
|---|---|---|
| `NMIQTREE.1.sh` | SLURM (4 nodes × 32 CPUs) | Maximum likelihood tree inference with ModelFinder (LG/WAG), ultrafast bootstrap (×1000), aLRT, and aBayes; no partition model |
| `NMIQTREE_mix5.sh` / `NMIQTREE_mix6.sh` | SLURM | IQ-TREE runs using mixture models (C60+LG+PMSF) with partition schemes |
| `AllBackbone.bp.merged.constraint.tre` | — | Backbone constraint tree used for constrained phylogenetic inference |

#### PhyloBayes (`pb_mpi/`)

| Script | Scheduler | Function |
|---|---|---|
| `NMBBtree1.PB.sh` / `NMBBtree4.PB.sh` | SLURM | PhyloBayes MPI inference (CAT-GTR) on backbone datasets of 58 or 43 taxa |
| `NMPB1.x.sh` / `NMPB4.x.sh` | SLURM | Worker scripts launching individual PhyloBayes MPI chains |
| `NMPBcomp.sh` | — | Convergence diagnostics: `bpcomp` and `tracecomp` run across all pairwise, triple, and quadruple chain combinations |
| `NMPBModelComparison.sh` | — | Cross-run model comparison for PhyloBayes analyses |

---

### Stage 6 — Molecular Clock (`6.MolecularClock/`)

This is the most computationally intensive stage. We ran MCMCtree using the approximate likelihood method across a large grid of analytical choices (two tree topologies, four calibration strategies, two clock models, two partition schemes, three rate priors) to explore how robust the divergence time estimates are. The R scripts then cluster the resulting trees and identify which analytical decisions most strongly influence the estimates.

#### Shell scripts

| Script | Scheduler | Function |
|---|---|---|
| `NMMCMCTr.8.sh` | — | Dispatcher: loops over all topology × calibration × clock × partition × rate combinations and submits MCMCtree jobs |
| `NMMCMCTr.8.x.1.sh` / `NMMCMCTr.8.x.2.sh` | SLURM | Workers: run the MCMCtree two-step approximation (Hessian/gradient computation, then MCMC posterior and prior sampling) |
| `NMMCMCTreeChecking.sh` | — | Checks chain convergence for all MCMCtree runs |
| `NMMCMCTreeSummary.8.sh` | — | Summarizes posterior time trees from completed MCMCtree runs |
| `NMMCMCTreeUncompressed.sh` | — | Locates the uncompressed MCMCtree output archives and re-compress them using a Zip compressor |
| `runPy.sh` | — | SBATCH wrapper for `SumMCMCtree.py` |

#### Python scripts

| Script | Function |
|---|---|
| `PAMLphylip.py` | Converts gene alignments and IQ-TREE partition files into PAML relaxed-PHYLIP format; supports multiple partition aggregation strategies (model-based, DBSCAN/KMeans/KMedoids clustering, tree-based) |
| `SumMCMCtree.py` | Reads MCMCtree MCMC samples and builds an ArviZ inference object; compares prior vs posterior node age distributions; outputs dated Nexus trees, summary tables, and convergence diagnostic plots |

#### R scripts

| Script | Function |
|---|---|
| `CalibrationPriorDistribution_6.r` | Translates fossil age constraints into user-specified statistical prior distributions (truncated Cauchy, skew-t, uniform); computes 95% HDI intervals; exports SVG plots and CSV tables for four calibration strategies |
| `TimeTreeAnalysesSelectionIntegration.r` | Integrates all MCMCtree results; parses convergence and model metadata from file names; uses a decision tree (RPART) to classify which analytical choices predict convergence |
| `TimeTreeGroupingEstimation.r` | Clusters converged time trees using PCA-DBSCAN, tSNE-DBSCAN, PRD-SBSCAN (based on custom-defiend distance), and autoencoder-DBSCAN to identify groups of similar time estimates and a consensus best estimate |
| `SensAna.r` | Visualizes sensitivity of key clade age estimates across all calibration strategies; produces streamer plots and annotated individual tree figures |
| `PlotSaveGroupTrees.r` | Builds averaged and centroid summary trees per cluster; computes 95% HPD intervals across trees; renders publication SVG figures with geological time scale overlays |
| `VariableAssessmentRDA.r` | Redundancy analysis (RDA) and variance partitioning to quantify how much of the variation in time estimates is explained by calibration strategy, clock model, substitution model, and partition scheme |
| `VariableAssessmentXGB.r` | XGBoost classification of convergence outcomes with SHAP importance values; ROC analysis; pinpoints the analytical parameters most responsible for non-convergence |

**Key R dependencies:** `ape`, `ggtree`, `treeio`, `vegan`, `ggplot2`, `xgboost`, `tidymodels`, `sn`, `tfprobability`, `doParallel`, `dbscan`, `Rtsne`, `torch`, `pracma`, `ggalt`

---

### Stage 7 — Ancestral State Reconstruction (`7.AncestralStateReconstruction/`)

The final analytical stage asks: when and how many times did parasitism evolve in nematodes? We fit HiSSE (Hidden State Speciation and Extinction) models to four independently dated phylogenies to account for uncertainty in both the tree topology and the divergence times. `SSEStateRateTest.r` then tests whether diversification rates differ between free-living and parasitic lineages through time.

| Script | Function |
|---|---|
| `HiSSEModels.r` | Core modelling function: fits the full HiSSE model suite (null, BiSSE-like, 2-state HiSSE, 4-state HiSSE, MiSSE); prunes redundant models; performs AIC/AICc selection; returns model-averaged rates and marginal ancestral state probabilities |
| `ACE.r` | Applies `HiSSEModels` across four dated phylogenies (cLGPMSF, C60PMSF, LG, traditional calibration); reconstructs ancestral free-living vs parasitic states; produces annotated tree figures; calls `SSEStateRateTest` for each tree |
| `SSEStateRateTest.r` | Bins phylogenetic nodes into geological time windows; categorizes nodes as FreeLiving / Parasitic / Hard / Soft; tests rate differences between states and before/after a temporal cut-point using ANCOVA and post-hoc multiple comparisons (please read the **in-line comments** for details about comparison methods!); generates rate-through-time plots and rate distribution histograms |

**Key R dependencies:** `hisse`, `ape`, `ggtree`, `treeio`, `geiger`, `ggplot2`, `patchwork`, `agricolae`, `car`, `plyr`, `dplyr`, `future`, `doParallel`, `readxl`, `geiger`, `phytools`

---

## Demo

The two demos below can be run on a standard desktop without any cluster access. They use lightweight inputs and complete in minutes, so they are a good way to check that your R/Python environments are set up correctly before committing to the full pipeline. Full data for the original study are available on Zenodo (10.5281/zenodo.19571088).

### Demo 1: Loci Selection

This filters a set of BUSCO alignments down to the best-quality loci for phylogenetics. You will need a `alignment_stats.csv` file produced by `NMAlignTrim.sh` (Stage 4), or you can use your own alignment quality table in the same format.

```bash
conda activate msaenv

python 4.MSA/LociSelection_aliba.py \
    -i msa/matrix/alignment_stats.csv \
    -o msa/selected_loci.txt
```

**Expected output:** A plain-text file with one file path per line, pointing to the loci that passed all quality thresholds — typically 400–600 out of ~800 BUSCO loci for a 100+ taxon dataset.

**Expected run time on a standard desktop:** Under 1 minute.

### Demo 2: Calibration Prior Visualization

This generates the fossil calibration prior distributions used in the molecular clock analysis. It is self-contained — just replace the ages according to your own time priors and run the script from its own directory.

```r
# In R (set working directory to 6.MolecularClock/ first):
setwd("6.MolecularClock")
source("calibration_prior_distribution_6.r")
```

**Expected output:** A set of SVG figures showing the prior probability density for each fossil calibration point, plus CSV tables with the distribution parameters and 95% HDI bounds.

**Expected run time on a standard desktop:** 5–15 minutes (the script uses `doParallel`; more cores will help).

---

## Instructions for Use

### Running the Full Pipeline on Your Data

If you are starting from scratch with a new taxon group, work through the steps below in order. Each step feeds the next. Where a step generates many parallel jobs (e.g., one per genome), this is handled automatically by the dispatcher scripts.

#### Step 1 — Collect sequence metadata from NCBI

Run the two Python scripts to search NCBI and build your download lists. Replace `"Nematoda"` with your target taxon:

```bash
python 1.SequenceAcquisition/get_assembly_info_20210819.py --input "Nematoda"
python 1.SequenceAcquisition/get_sra_info_max20210831.py --input "Nematoda"
```

Both scripts output Excel files. Review and curate these manually to remove assemblies with very low coverage, redundant strains, or unsuitable library types before proceeding.

#### Step 2 — Download the data

```bash
sbatch 1.SequenceAcquisition/GenomeDownload.sh assembly_info.csv
```

The script runs up to 20 parallel downloads and writes its progress to `CLine.log`, so if it is interrupted you can resubmit and it will pick up where it left off. The CSV file `Maker-EST.csv` in the same directory shows an example of the `assembly_info.csv` in the above command.

#### Step 3 — Assemble transcriptomes (RNA-seq)

Submit one job per sample. The script handles quality control, Trinity assembly, and TransDecoder ORF prediction in a single run:

```bash
sbatch 2.Assembly/NMTrinity.9.x.sh <sample_name>
```

If your cluster supports submitting multiple jobs simultaneously, it is recommended to use `NMTrinity.9.sh` rather than executing `NMTrinity.9.x.sh` scripts sequentially.

#### Step 4 — Assemble genomes (WGS)

```bash
qsub 2.Assembly/Nematoda_spades.0.pbs
```

To run these jobs in parallel, you can enable parallel processing by uncommenting the specified lines. However, doing this on a single node is not recommended. A better approach is to submit multiple *Slurm* jobs, running each on a separate node.

#### Step 5 — Decontaminate assemblies

First screen for vector contamination, then remove non-target taxonomic sequences. Finally, use the pie chart script to visually confirm the taxonomic composition of what remains:

```bash
# Remove vector-contaminated regions
python 2.Assembly/RmVecRNA.mp.py --fasta assembly.fasta --blast univec_hits.csv

# Remove non-nematode sequences (rank-based)
python 2.Assembly/DeContaminRNA-rk.py --fasta assembly.fasta --blast diamond_hits.csv --valid Nematoda

# Visualize the remaining taxonomic composition
python 2.Assembly/SeqOriginPiePlot.py --fasta assembly.fasta --blast diamond_hits.csv
```

The preferred method for executing these scripts is to run the following command:

```bash
sbatch 2.Assembly/UniVec.Decontamination.6.sh
```

#### Step 6 — Gene prediction

Repeat masking should be performed prior to running MAKER, as MAKER utilizes the masked genome by default. In addition to MAKER, MetaEuk is also employed for gene prediction. The `NMPredMaker.6.a.x.sh` script integrates Repeat Masking, MAKER, and MetaEuk, and can be executed for each genomic assembly sequentially or as part of multiple jobs submitted simultaneously using `NMPredMaker.6.a.sh`:

```bash
# Run MAKER gene prediction (dispatches one SBATCH job per genome automatically)
sbatch 3.Prediction/NMPredMaker.6.a.sh
```

#### Step 7 — Consolidate BUSCO orthologs across all taxa

Copy or move the predicted amino acid sequences (*.pep.faa) to the directory `$DIR_current/pep`, and then submit a Slurm job using the following command:

```bash
sbatch 3.Prediction/NMBusco.sh
```

This process gathers the BUSCO protein sequences for each species and generates a single unified FASTA file for each locus, preparing them for alignment. Typically, a single job on one node with approximately 32 vCPUs is sufficient, and it runs quite quickly. If you have a large volume of data, you can utilize as many nodes as permitted, but running in a single job is proper than distributing several jobs.

#### Step 8 — Multiple sequence alignment and loci selection

```bash
# Align, trim, and score all BUSCO loci in parallel
sbatch 4.MSA/NMAlignTrim.sh

# Filter loci by quality (quick approach) and using gene-tree evidence too, then concatenate selected loci into a single supermatrix
sbatch 4.MSA/NMFilterConcat.sh
```

#### Step 9 — Phylogenetic inference

We recommend running both IQ-TREE and PhyloBayes and comparing the topologies before proceeding to dating:

```bash
# Maximum likelihood tree, it results a best fitting substitution model and an initial ML tree that can be used to identify fully supported clades for subsequent clade-wise subsampling
sbatch 5.Phylogenetics/IQTREE/NMIQTREE.1.sh

# Maximum likelihood trees with more complicated substitution models, LG4X and Poisson-C60 and LG-C60-PMSF
sbatch 5.Phylogenetics/IQTREE/NMIQTREE_mix5.sh

# Bayesian inference (runs multiple parallel chains)
sbatch 5.Phylogenetics/pb_mpi/FBB_58taxa/NMBBtree1.PB.sh

# Bayesian inference (runs multiple parallel chains)
sbatch 5.Phylogenetics/pb_mpi/FBB_43taxa/NMBBtree4.PB.sh

# Verify if the chains have converged (this should be done after adequate sampling). As indicated in the example command, check chains 1 to 4 while discarding the first 1000 generations as burn-in.
sbatch 5.Phylogenetics/pb_mpi/NMPBcomp.sh 1 4 1000

# Maximum likelihood tree using PMSF which was guided by a LG tree that was caonstrained by a combined Bayesian tree from the two Bayesian phylogenies above.
sbatch 5.Phylogenetics/IQTREE/NMIQTREE_mix6.sh
```

#### Step 10 — Molecular clock dating

This step requires some manual preparation: edit `mcmctree.ctl` to specify your tree, alignment, and calibration nodes before running the dispatcher. The input sequence file can be the 184-taxa output (full dataset with outgroups) of `4.MSA\NMFilterConcat.sh`.

```bash
# Convert alignment to PAML relaxed-PHYLIP format
python 6.MolecularClock/PAMLphylip.py --seq FullData184taxa171097AA.phy --partitions partitions.nex --outdir timing/

# Visualize and set fossil calibration priors (R), replace the ages with your time priors before use
Rscript 6.MolecularClock/calibration_prior_distribution_6.r

# Submit all MCMCtree runs with a fixed sample fraction(loops over all analytical combinations), if you want to test different sample fractions, copy the whole directory and modify the `mcmctree.ctl` file and run the following command again:
sbatch 6.MolecularClock/NMMCMCTr.8.sh

# Summarize the posteriors and priors once the runs are complete. You don't need to wait for all jobs to finish, as the script is designed to identify which prior and posterior MCMC sampling chains have all completed first.
sbatch 6.MolecularClock/NMMCMCTreeSummary.8.sh

# Cluster the resulting time trees and identify robust estimates
Rscript 6.MolecularClock/TimeTreeAnalysesSelectionIntegration.r
Rscript 6.MolecularClock/TimeTreeGroupingEstimation.r

# Assess which analytical choices most influence the results
Rscript 6.MolecularClock/VariableAssessmentRDA.r
Rscript 6.MolecularClock/VariableAssessmentXGB.r
```

Conduct sensitivity analysis across calibration strategies. Before this, you must first perform molecular clock analyses (`NMMCMCTr.8.sh`) using sensitivity data and under the specified sensitivity conditions.

```bash
Rscript 6.MolecularClock/SensAna.r
```

#### Step 11 — Ancestral state reconstruction

Ensure that all libraries and dependencies are installed and available before executing the script `ACE.r`. This script will load `HiSSEModels.r` and `SSEStateRateTest.r` prior to processing.

```bash
Rscript 7.AncestralStateReconstruction/ACE.r
```

---

## Expected Outputs

| Stage | Key Output Files |
|---|---|
| 1. Acquisition | `assembly_info.xlsx`, `sra_info.xlsx`, downloaded FASTA/SRA files |
| 2. Assembly | `*_trinity.Trinity.fasta`, `*_scaffolds.fasta`, `*.tdc.pep.faa` |
| 3. Prediction | `*.final.maker.pep.faa`, `*.mtek.ur90.pep.faa`, BUSCO summary reports |
| 4. MSA | `*.mafft.fasta`, `*.bmge.fasta`, `alignment_stats.csv`, `selected_loci.txt`, `supermatrix.phy` |
| 5. Phylogenetics | `*.contree` (ML tree), `*.treefile` (ML tree),`*.pb` chain files, `bpcomp*.log` (convergence) |
| 6. Molecular Clock | `*_posterior_tree_good.tre`, `*_posterior_tree_acpt.tre`, `*_posterior_tree_bad.tre`, `KeyCladeTime_SATrees.csv`, calibration prior SVG plots, `AllPosteriorTimeTree.*.RData` |
| 7. Ancestral State | `HiSSE_ModelSelection-*.csv`, `HiSSE_RateData*-*.csv`, `ANCOVA report*.txt`, `Rate Change Through Time-*.svg`, `Rate Distributions Histogram-*.svg` |

---

## Dependencies Summary

### Python Packages

| Package | Version (tested) | Used In |
|---|---|---|
| `biopython` | ≥ 1.79 | Stages 1, 2, 3, 4, 6 |
| `pandas` | ≥ 1.3 | Stages 1, 2, 3, 4, 6 |
| `numpy` | ≥ 1.21 | Stages 2, 3, 4, 6 |
| `matplotlib` | ≥ 3.4 | Stages 2, 3, 4 |
| `beautifulsoup4` | ≥ 4.9 | Stage 1 |
| `tqdm` | ≥ 4.60 | Stage 1 |
| `scikit-learn` | ≥ 1.0 | Stage 6 |
| `sklearn-extra` | ≥ 0.2 | Stage 6 |
| `kneed` | ≥ 0.7 | Stage 6 |
| `arviz` | ≥ 0.12 | Stage 6 |
| `seaborn` | ≥ 0.11 | Stage 6 |
| `xarray` | ≥ 0.19 | Stage 6 |
| `scipy` | ≥ 1.7 | Stage 6 |
| `argparse` | stdlib | Throughout |
| `multiprocessing` | stdlib | Stage 2 |

### R Packages

| Package | Used In |
|---|---|
| `ape` | Stages 6, 7 |
| `ggtree` / `treeio` | Stages 6, 7 |
| `ggplot2` / `patchwork` | Stages 6, 7 |
| `hisse` | Stage 7 |
| `geiger` | Stage 7 |
| `vegan` | Stage 6 |
| `xgboost` | Stage 6 |
| `tidymodels` / `recipes` / `parsnip` | Stage 6 |
| `sn` | Stage 6 |
| `tfprobability` | Stage 6 |
| `doParallel` / `future` | Stages 6, 7 |
| `dbscan` / `Rtsne` / `fpc` | Stage 6 |
| `torch` / `brulee` | Stage 6 |
| `rpart` / `caret` | Stage 6 |
| `car` / `agricolae` | Stage 7 |
| `phytools` / `ips` | Stage 7 |
| `readxl` | Stage 7 |
| `pracma` | Stage 6 |
| `ggalt` / `grid` | Stage 6 |

### Bioinformatics Tools (installed via Conda)

| Tool | Version | Stage |
|---|---|---|
| SPAdes | 3.15.0 | 2 |
| SOAPdenovo2 / SOAPdenovo-Trans | latest | 2 |
| Trinity | 2.13–2.15 | 2 |
| TransDecoder | latest | 2 |
| AfterQC | latest | 2 |
| BBMap / BBtools | latest | 2 |
| RepeatMasker / RepeatModeler | 4.2.3 | 3 |
| MAKER | latest | 3 |
| MetaEuk | 7.x | 3 |
| BUSCO | 6.0.0 | 3 |
| BLAST+ | 2.17.0 | 2, 3 |
| DIAMOND | 2.1.24 | 2, 3 |
| MAFFT | latest | 4 |
| BMGE | 1.12 | 4 |
| PhyKit | latest | 4 |
| IQ-TREE2 | 2.1.3 | 5 |
| PhyloBayes MPI | latest | 5 |
| PAML / MCMCtree | latest | 6 |

---

## Notes and Common Gotchas

- **Absolute paths** — Before running anything, do a global search-and-replace of `/To/Your/Directory/` with your own base directory. A quick `grep -r 'To/Your/Directory'` from the repo root will confirm all instances.
- **SLURM / PBS parameters** — job resource requests (nodes, CPUs, memory, partition names) are tuned for our specific cluster. Adjust them to fit your own scheduler and queue policies.
- **NCBI email** — the Python scripts that query NCBI require a registered email address. Replace `A.N.Other@example.com` or `lianges.luex@gmail.com (my address)` with your own wherever it appears.
- **`Classification.csv`** — this file (in `6.MolecularClock/`) maps taxon names to their phylogenetic and taxonomic groups. It is required by several scripts in Stages 6 and 7. If you are using a different taxon set, you will need to update this file accordingly.
- **`mcmctree.ctl`** — the PAML control file needs to be manually edited to point to your alignment, your tree file, and your calibration nodes before launching any MCMCtree jobs. The current file reflects our specific analysis and will not work out-of-the-box for a different dataset.
- **R package `hisse`** — install from GitHub for the most recent version: `devtools::install_github("thej022214/hisse")`.
