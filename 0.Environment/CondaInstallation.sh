#!/bin/sh

# Update the Conda package manager and the Anaconda distribution.
conda update conda
conda update anaconda

# Configure standard bioinformatics channels. Order matters for package resolution.
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels thiesgehrmann

# Install core Python and scientific packages in the base environment.
# Use Python 3.10 for compatibility with current Conda packages.
conda install python=3.10
conda install biopython
conda install arviz=1.0.0 # need testing for compatibility with SumMCMCtree.py
conda install numba
conda install -c conda-forge scikit-learn-extra kneed

# Clean unused packages and cache to free disk space.
conda clean -p && conda clean -y --all

# List available Conda environments in example installation paths.
conda env list /To/Your/Directory/anaconda3
conda env list /public1/home/scb4616/apps/anaconda3

# Create a BLAST environment for sequence search tools.
conda create -n blastenv
conda install -n blastenv blast=2.17.0 diamond=2.1.24

# Create a MAKER annotation environment with gene prediction and evidence integration tools.
conda create -n makerenv
conda install -n makerenv maker
conda install -c thiesgehrmann -n makerenv genemark_es
conda install -n makerenv star bedtools repeatmodeler pasa metaeuk

# Create a repeat-masking environment for transposable element annotation and repeat libraries.
conda create -n repmaskenv python=3.7
conda install -n repmaskenv repeatmasker=4.2.3 repeatmodeler edta dfam

# Create a BRAKER2 environment for automated gene prediction using RNA-seq evidence.
conda create -n braker2env
conda install -n braker2env braker2

# Create a BUSCO environment for genome completeness assessment and plotting support.
# Newer BUSCO versions require Python 3 and are best run with a modern Python 3.10 environment.
conda create -n buscoenv python=3.10
conda install -n buscoenv r-base r-ggplot2 metaeuk=7.bba0d80 prodigal sepp busco=6.0.0

# Create an OrthoFinder environment for orthology detection, plus UCX support for fast transfers.
# OrthoFinder v3.x is Python 3 compatible; Python 2 is only needed for legacy OrthoFinder 2.x workflows.
conda create -n orthofinderenv python=3.10
conda install -n orthofinderenv orthofinder
conda install -n orthofinderenv -c conda-forge ucx

# Create a genome browser and visualization environment.
conda create -n gmshowenv
conda install -n gmshowenv jbrowse2 apollo pygenometracks igv

# Create a Python 3.12 environment for the newest compatible packages.
conda create -n newpy python=3.12

# Create an assembly environment for genome assembly tools.
conda create -n assemblingenv
conda install -n assemblingenv spades soapdenovo2 soapdenovo2-gapcloser soapdenovo2-prepare

# Create an RNA-seq environment for transcriptome assembly.
conda create -n rnaseqenv
conda install -n rnaseqenv trinity=2.15.2 soapdenovo-trans transdecoder
# soapdenovo2 can also be used for transcriptome assembly, but it may have some issues.

# Example TransRate environment setup for transcriptome assembly quality evaluation (disabled by default).
# conda create -n transrateenv
# conda install -n transrateenv transrate transrate-tools

# Create a sequence quality control environment.
conda create -n seqqc
conda install -n seqqc fastp fastqc besst

# Create a multiple sequence alignment environment.
conda create -n msaenv
conda install -n msaenv mafft bmge
conda install -n msaenv -c jlsteenwyk phykit clipkit

# Create a PhyloBayes environment for Bayesian phylogenetic inference.
conda create -n phylobayes
conda install -n phylobayes -c bioconda phylobayes-mpi=1.9

# Show current Conda environments and environment info.
conda env list
conda info -e
# Example commands for activating and removing Conda environments.
source activate your_env_name
conda remove -n your_env_name --all
conda remove --name $your_env_name $package_names