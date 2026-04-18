#!/To/Your/Directory/anaconda3/bin/python3.9

"""
LociSelection_aliba.py - Select High-Quality Loci for Phylogenetic Analysis

This script filters loci (gene alignments) based on multiple quality metrics to select
the most informative sequences for downstream phylogenetic analysis. It uses quantile-based
thresholds to ensure selected loci meet minimum quality standards.

The script evaluates loci on several criteria:
- ngen: Number of genes/sequences in the alignment
- pis: Percentage of parsimony-informative sites
- tot: Total alignment length
- RCV: Relative Composition Variability (measure of compositional heterogeneity)
- erate: Error/missing data rate

Loci must meet or exceed the 50th percentile for ngen, pis, and tot, while staying
below the 90th percentile for RCV and erate.

Input: CSV file with alignment statistics (output from NMAlignTrim.sh)
Output: Text file with paths to selected loci, one per line

Dependencies: numpy, pandas
"""

import numpy as np
import pandas as pd
import argparse

# Set up command-line argument parser
parser = argparse.ArgumentParser(
    description="Filter loci based on quality metrics for phylogenetic analysis",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Input CSV format (comma-separated):
fastapath,ngen,pis,tot,ppis,GC,RCV,erate
/path/to/locus.fasta,25,0.15,1200,0.12,0.45,0.02,0.05

Selection criteria:
- ngen >= 50th percentile (sequence count)
- pis >= 50th percentile (parsimony-informative sites)
- tot >= 50th percentile (alignment length)
- RCV <= 90th percentile (composition variability)
- erate <= 90th percentile (error rate)

Example:
  python LociSelection_aliba.py -i alignment_stats.csv -o selected_loci.txt
    """
)
parser.add_argument('--infile', '-i',
                    type=str,
                    required=True,
                    help='Input CSV file with alignment statistics (output from NMAlignTrim.sh)')
parser.add_argument('--outfile', '-o',
                    type=str,
                    required=True,
                    help='Output text file containing paths to selected loci')

args = parser.parse_args()

# Process file paths (handle Windows backslashes)
infile = args.infile.replace('\\','/')
outfile = args.outfile.replace('\\','/')

# Read alignment statistics from CSV file
# Columns: fastapath, ngen, pis, tot, ppis, GC, RCV, erate
LociFilterMeasure = pd.read_table(infile, index_col=None, 
                                 names=['fastapath', 'ngen', 'pis', 'tot', 'ppis', 'GC', 'RCV', 'erate'], 
                                 sep=",", engine='python')

# Apply quality filters using quantile-based thresholds
# Select loci that meet all criteria:
# - At least median number of sequences (ngen >= 50th percentile)
# - At least median parsimony-informative sites (pis >= 50th percentile)  
# - At least median alignment length (tot >= 50th percentile)
# - Below 90th percentile for composition variability (RCV <= 90th percentile)
# - Below 90th percentile for error rate (erate <= 90th percentile)
SelectedLociPath = LociFilterMeasure.loc[
    (LociFilterMeasure['ngen'] >= np.quantile(LociFilterMeasure['ngen'], 0.5)) &  # Sequence count filter
    (LociFilterMeasure['pis'] >= np.quantile(LociFilterMeasure['pis'], 0.5)) &    # Parsimony-informative sites filter
    (LociFilterMeasure['tot'] >= np.quantile(LociFilterMeasure['tot'], 0.5)) &    # Alignment length filter
    (LociFilterMeasure['RCV'] <= np.quantile(LociFilterMeasure['RCV'], 0.9)) &    # Composition variability filter
    (LociFilterMeasure['erate'] <= np.quantile(LociFilterMeasure['erate'], 0.9)),  # Error rate filter
    'fastapath'  # Extract only the file paths
]

# Write selected locus paths to output file (one path per line, no header)
SelectedLociPath.to_csv(outfile, index=False, header=False)

print(f"Selected {len(SelectedLociPath)} loci out of {len(LociFilterMeasure)} total loci based on quality filters.")

