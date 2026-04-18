#!/bin/sh

# Initialization script for the Nematoda analysis workspace.
# Set the path to the Python interpreter used for downstream scripts.
PythonIntepreter=

# Root directory for the Anaconda installation if needed by scripts.
AnacondaDir=

# Default directory for downloaded or generated sequence databases.
DatabaseDIR=

# Base directory for the project workspace. This should be set before running init.
BaseinDir=

# Create the project directory hierarchy for data, analysis, and results.
mkdir -p "$BaseinDir"
mkdir -p "$BaseinDir/data"
mkdir -p "$BaseinDir/data/assembly"
mkdir -p "$BaseinDir/data/rnaseq"
mkdir -p "$BaseinDir/data/wgs"
mkdir -p "$BaseinDir/data/preds"
mkdir -p "$BaseinDir/data/preds/augustus_config"
mkdir -p "$BaseinDir/data/busco"
mkdir -p "$BaseinDir/msa"
mkdir -p "$BaseinDir/msa/align"
mkdir -p "$BaseinDir/msa/trim"
mkdir -p "$BaseinDir/msa/matrix"
mkdir -p "$BaseinDir/phylo"
mkdir -p "$BaseinDir/timing"
