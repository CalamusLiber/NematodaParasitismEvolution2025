#!/To/Your/Directory/anaconda3/bin/python3.9

"""
LociSelection_treba.py - Advanced BUSCO Loci Selection and Concatenation

This script performs sophisticated filtering and concatenation of BUSCO gene alignments
for phylogenetic analysis. It integrates multiple quality metrics from both sequence
alignments and gene trees to select the most reliable loci.

The script evaluates loci based on:
Alignment metrics: nGene, nPIS, TrimAlignLength, pPIS, GC, RCV, EvoRate
Gene tree metrics: AvgBoots, Treeness, RCV, TOR, DVMC

Optional spurious homolog detection removes potentially problematic sequences
identified as outliers in gene trees.

Workflow:
1. Filter loci based on quality thresholds
2. Remove spurious homologs (if provided)
3. Concatenate selected alignments into supermatrix
4. Generate statistics and visualizations

Input: Quality metric files from NMAlignTrim.sh and NMFilterConcat.sh
Output: Concatenated alignments in PHYLIP/FASTA format, partition files, statistics

Dependencies: Biopython, matplotlib, pandas, numpy
"""

import numpy as np
import pandas as pd
import os, re
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt
from matplotlib import axis, gridspec as gs
import argparse

def DropSeq(alignpath, seqstodrop, exclude, outdir, mingenes=1):
    """
    Remove spurious sequences from an alignment and write filtered version.
    
    Removes sequences identified as spurious homologs from gene tree analysis,
    along with any user-specified taxa to exclude.
    
    Args:
        alignpath (str): Path to input alignment file
        seqstodrop (pd.DataFrame): Sequences to remove (from spurious homolog detection)
        exclude (pd.Series): Additional taxa to exclude
        outdir (str): Output directory
        mingenes (int): Minimum number of sequences required to keep alignment
        
    Returns:
        pd.Series: Sequences in the filtered alignment, or empty series if too few sequences
    """
    alignname = alignpath.split('/')[-1].split('.')[0]
    
    # Determine output directory
    if outdir == '':
        outdir = alignpath.split(alignname)[0][:-1]
    
    # Read alignment
    Align = AlignIO.read(alignpath, 'fasta')
    
    # Combine sequences to drop (spurious + excluded)
    seqstodrop_combined = pd.concat([pd.Series(seqstodrop.loc[alignname, 'species']), exclude]) if alignname in seqstodrop.index else exclude
    
    # Filter out unwanted sequences
    Align = [seq for seq in Align if seq.id not in set(seqstodrop_combined.values)]
    
    # Write filtered alignment if it meets minimum requirements
    if len(Align) >= mingenes:
        Align = MultipleSeqAlignment(Align)
        AlignIO.write(Align, outdir + '/' + alignname + '.trbdp.rmsp.faa', 'fasta')
        return pd.Series([seq.seq for seq in Align], index=[seq.id for seq in Align], name=alignname)
    else:
        return pd.Series([], dtype=str)

def ReadSeq(alignpath, exclude, outdir, mingenes=1):
    """
    Read and filter alignment, removing excluded taxa.
    
    Similar to DropSeq but without spurious homolog removal.
    
    Args:
        alignpath (str): Path to input alignment file
        exclude (pd.Series): Taxa to exclude from analysis
        outdir (str): Output directory
        mingenes (int): Minimum number of sequences required
        
    Returns:
        pd.Series: Sequences in the filtered alignment
    """
    alignname = alignpath.split('/')[-1].split('.')[0]
    
    # Read alignment
    Align = AlignIO.read(alignpath, 'fasta')
    
    # Filter out excluded taxa
    Align = [seq for seq in Align if seq.id not in exclude.values]
    
    # Write filtered alignment if it meets requirements
    if len(Align) >= mingenes:
        Align = MultipleSeqAlignment(Align)
        AlignIO.write(Align, outdir + '/' + alignname + '.trbdp.faa', 'fasta')
        return pd.Series([seq.seq for seq in Align], index=[seq.id for seq in Align], name=alignname)
    else:
        return pd.Series([], dtype=str)

def ConcatFas2Phy(alignseries, outfile, partitionfile=True, returnaspyobj=True):
    """
    Concatenate multiple alignments into a supermatrix and generate partition file.
    
    Takes a series of alignments, sorts them by length (longest first), fills missing
    taxa with gaps, and concatenates into a single supermatrix. Generates PHYLIP,
    FASTA, and partition files for downstream analysis.
    
    Args:
        alignseries (pd.Series): Series of alignments to concatenate
        outfile (str): Base output filename (without extension)
        partitionfile (bool): Whether to generate partition file
        returnaspyobj (bool): Whether to return alignment statistics
        
    Returns:
        pd.DataFrame: Alignment integrity statistics (if returnaspyobj=True)
    """
    # Sort alignments by length (longest first)
    lenlist = alignseries.apply(lambda x: len(x[0]))
    lenlist.sort_values(inplace=True, ascending=False)
    
    # Create missing sequence placeholders
    missingseq = lenlist.apply(lambda x: Seq('?' * x))
    
    # Reorder alignments by length
    alignseries = alignseries[lenlist.index]
    
    # Concatenate alignments into dataframe
    alignframe = pd.concat(alignseries.tolist(), axis=1, join='outer', ignore_index=False)
    
    # Fill missing taxa with gaps
    for x in alignframe.columns:
        for y in alignframe.loc[alignframe.loc[:, x].isna(), x].index:
            alignframe.loc[y, x] = missingseq[x]
    
    # Sort taxa alphabetically
    alignframe.sort_index(axis=0, inplace=True)
    
    # Create Biopython alignment object
    alignments = MultipleSeqAlignment(
        alignframe.apply(lambda x: SeqRecord(Seq(''.join([str(s) for s in x])), id=x.name), axis=1).tolist()
    )
    
    # Write output files
    AlignIO.write(alignments, outfile + '.phy', 'phylip-relaxed')
    AlignIO.write(alignments, outfile + '.faa', 'fasta')
    
    # Generate partition file for partitioned analysis
    if partitionfile:
        xlength = pd.concat([pd.Series(0), lenlist]).cumsum()
        ptflines = ['  charset ' + xlength.index[i+1] + ' = ' + str(xlength.iloc[i]+1) + '-' + str(xlength.iloc[i+1]) + ';\n' 
                   for i in range(len(xlength)-1)]
        with open(outfile + '.partition.nex', mode='w', newline='\n') as ptf:
            ptf.write('#nexus\nbegin sets;\n')
            ptf.writelines(ptflines)
            ptf.write('end;\n')
    
    # Calculate and return alignment statistics
    if returnaspyobj:
        totlen = sum(lenlist)
        totspp = alignframe.shape[0]
        
        # Calculate integrity statistics
        alignlength = alignframe.applymap(lambda x: sum([p.isalnum() for p in x]))
        alignlengthout = alignlength.copy()
        
        # Per-taxon statistics
        alignlengthout.loc[:, 'Integrity(%,total=' + str(totlen) + ')'] = alignlength.apply(
            lambda x: round(sum(x) / totlen * 100, 2), axis=1)
        alignlengthout.loc[:, 'NumGenes'] = alignlength.apply(lambda x: sum([i > 0 for i in x]), axis=1)
        
        # Per-gene statistics
        alignlengthout.loc['Integrity(%)', :] = round(
            alignlength.apply(lambda x: sum(x)/totspp, axis=0) / lenlist * 100, 2)
        alignlengthout.loc['MissingGenes', :] = alignlength.apply(lambda x: sum([i == 0 for i in x]), axis=0)
        
        return alignlengthout

def HeatMap(taxastat, outfile, dpi=300, nfsize=50, tickstep=10):
    """
    Generate a three-panel heatmap visualization of alignment integrity.
    
    Creates a figure with three subplots showing:
    1. Log-transformed ungapped sequence lengths per gene per taxon
    2. Binary presence/absence of genes per taxon
    3. Bar chart of missing genes per locus
    
    Args:
        taxastat (pd.DataFrame): Alignment statistics matrix
        outfile (str): Output filename for the plot
        dpi (int): Resolution for saved image
        nfsize (int): Figure size multiplier
        tickstep (int): Step size for axis tick labels
    """
    # Prepare data for plotting
    ts_log = taxastat.iloc[:-2, :-2].apply(np.log).fillna(0)  # Log lengths, fill NA with 0
    ts_bin = taxastat.iloc[:-2, :-2].applymap(lambda x: x > 0)  # Binary presence/absence
    
    # Set up tick positions
    xlabpos = np.arange(len(ts_log.columns), step=tickstep*3)
    ylabpos = np.arange(len(ts_log.index), step=tickstep)
    
    # Create figure with three subplots
    fig = plt.figure(figsize=(3*nfsize, 2*nfsize))
    gsp = gs.GridSpec(3, 1, height_ratios=[1, 1, 1])
    
    # Top plot: Log ungapped lengths
    axes0 = plt.subplot(gsp[0])
    htp0 = axes0.imshow(ts_log, cmap=plt.cm.jet, aspect='auto')
    axes0.set_title('Log ungappy length of gene vs. taxon', fontdict={'fontsize':25})
    axes0.set_xticks([])
    axes0.set_yticks(ylabpos)
    axes0.set_yticklabels(labels=ts_log.index[ylabpos].tolist(), fontdict={'fontsize':3})
    axes0.set_ylabel('Taxon', fontdict={'fontsize':18})
    
    # Middle plot: Gene presence/absence
    axes1 = plt.subplot(gsp[1])
    htp1 = axes1.imshow(ts_bin, cmap=plt.cm.gray, aspect='auto')
    axes1.set_title('Presence of gene vs. taxon', fontdict={'fontsize':25})
    axes1.set_xticks([])
    axes1.set_yticks(ylabpos)
    axes1.set_yticklabels(labels=ts_bin.index[ylabpos].tolist(), fontdict={'fontsize':3})
    axes1.set_ylabel('Taxon', fontdict={'fontsize':18})
    
    # Bottom plot: Missing genes bar chart
    axes2 = plt.subplot(gsp[2])
    axes2.bar(taxastat.iloc[-1, 2:].index, taxastat.iloc[-1, 2:], width=1, color='orange')
    axes2.margins(x=0)
    axes2.set_xlabel('Trimmed Alignment of BUSCO Genes', fontdict={'fontsize':18})
    axes2.set_ylabel('Missing', fontdict={'fontsize':18})
    axes2.set_title('Count of missing genes', fontdict={'fontsize':25})
    axes2.set_xticks(xlabpos)
    axes2.set_xticklabels(labels=ts_bin.columns[xlabpos].tolist(), fontdict={'fontsize':3})
    plt.setp(axes2.get_xticklabels(), rotation=270)
    
    # Finalize and save plot
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    fig.savefig(outfile, dpi=dpi)

if __name__ == "__main__":
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Advanced BUSCO loci selection and concatenation for phylogenetics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input Files:
  -a: Alignment metrics CSV (from NMAlignTrim.sh)
       Columns: fastapath,nGene,nPIS,TrimAlignLength,pPIS,GC,RCV,EvoRate
  -g: Gene tree metrics CSV (from NMFilterConcat.sh)  
       Columns: genetree,AvgBoots,Treeness,RCV,TOR,DVMC
  -s: Spurious homologs CSV (optional, from NMFilterConcat.sh)
       Columns: genetree,species,termbranchlength,threshold,medianbranchlength

Selection Criteria:
  - AvgBoots >= 10th percentile
  - TOR >= 10th percentile  
  - DVMC <= 90th percentile

Example:
  python LociSelection_treba.py -a align_metrics.csv -g tree_metrics.csv -o results/
  python LociSelection_treba.py -a align.csv -g tree.csv -s spurious.csv -e "taxon1,taxon2" -o output/
        """
    )
    parser.add_argument('--alignmentmeasure', '-a',
                        type=str,
                        required=True,
                        help='CSV file with alignment quality metrics from NMAlignTrim.sh')
    parser.add_argument('--genetreemeasure', '-g',
                        type=str,
                        required=True,
                        help='CSV file with gene tree quality metrics from NMFilterConcat.sh')
    parser.add_argument('--spurioushomolog', '-s',
                        type=str,
                        default='',
                        help='Optional CSV file with spurious homolog detections')
    parser.add_argument('--outdir', '-o',
                        type=str,
                        default='./',
                        help='Output directory for results')
    parser.add_argument('--exclude', '-e',
                        type=str,
                        default='',
                        help='Comma-separated list of taxa to exclude from analysis')
    
    args = parser.parse_args()

    # Process input paths
    alignmentmeasure = args.alignmentmeasure.replace('\\','/')
    genetreemeasure = args.genetreemeasure.replace('\\','/')
    outdir = args.outdir.replace('\\','/')
    curdir = os.getcwd().replace('\\','/')
    
    # Convert to absolute paths
    if not os.path.isabs(outdir):
        outdir = curdir + '/' + outdir.replace('./','')
    
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir, mode=0o777, exist_ok=True)
    
    # Parse excluded taxa
    if args.exclude != '':
        excl = pd.Series(args.exclude.split(','), dtype=str)
    else:
        excl = pd.Series([], dtype=str)
    
    # Read alignment quality metrics
    AlignmentMeasure = pd.read_table(alignmentmeasure, index_col=None, 
                                   names=['fastapath', 'nGene', 'nPIS', 'TrimAlignLength', 'pPIS', 'GC', 'RCV', 'EvoRate'], 
                                   sep=",", engine='python')
    AlignmentMeasure.index = AlignmentMeasure['fastapath'].apply(lambda x: x.split('/')[-1].split('.')[0])
    
    # Read gene tree quality metrics
    GenetreeMeasure = pd.read_table(genetreemeasure, index_col=None, 
                                   names=['genetree', 'AvgBoots', 'Treeness', 'RCV', 'TOR', 'DVMC'], 
                                   sep=",", engine='python')
    GenetreeMeasure.index = GenetreeMeasure['genetree'].apply(lambda x: x.split('/')[-1].split('.')[0])
    
    # Combine alignment and tree metrics
    AdwiseLociMeasure = pd.concat([AlignmentMeasure, GenetreeMeasure], axis=1, join='inner', ignore_index=False)
    
    # Apply quality filters to select loci
    # Criteria: high bootstrap support, tree resolution, low deviation from molecular clock
    SelectedLoci = AdwiseLociMeasure.loc[
        (AdwiseLociMeasure['AvgBoots'] >= np.quantile(AdwiseLociMeasure['AvgBoots'], 0.1)) &  # Bootstrap support >= 10th percentile
        (AdwiseLociMeasure['TOR'] >= np.quantile(AdwiseLociMeasure['TOR'], 0.1)) &            # Tree resolution >= 10th percentile
        (AdwiseLociMeasure['DVMC'] <= np.quantile(AdwiseLociMeasure['DVMC'], 0.9)),            # Clock deviation <= 90th percentile
        :
    ].sort_index()
    
    # Export selected loci statistics
    SelectedLoci[['nGene', 'TrimAlignLength', 'nPIS', 'pPIS', 'RCV', 'EvoRate', 'AvgBoots', 'Treeness', 'TOR', 'DVMC']].to_excel(
        outdir + '/FinalSelectedGenes.xlsx')
    
    # Process alignments with or without spurious homolog removal
    if args.spurioushomolog != '':
        spurioushomolog = args.spurioushomolog.replace('\\','/')
        if os.path.exists(spurioushomolog):
            # Read spurious homolog data
            SpuriousHomolog = pd.read_table(spurioushomolog, index_col=None, 
                                          names=['genetree', 'species', 'termbranchlength', 'threshold', 'medianbranchlength'], 
                                          sep=",", engine='python')
            SpuriousHomolog.index = SpuriousHomolog['genetree'].apply(lambda x: x.split('/')[-1].split('.')[0])
            
            # Process alignments with spurious sequence removal
            dq = lambda x: DropSeq(x, SpuriousHomolog, exclude=excl, outdir=outdir, mingenes=1)
            AlignSeries = list(map(dq, SelectedLoci['fastapath']))
            AlignSeries = pd.Series(AlignSeries, index=[x.name for x in AlignSeries])
            AlignSeries = AlignSeries[AlignSeries.apply(len) >= 10]  # Require at least 10 taxa
            
            # Concatenate and analyze
            AlignLength = ConcatFas2Phy(AlignSeries, outfile=outdir + '/ConcatenatedMatrix.rmsp', 
                                       partitionfile=True, returnaspyobj=True)
            AlignLength.to_excel(outdir + '/MultipleSequenceAlignmentsLogLength.rmsp.xlsx')
            HeatMap(AlignLength, outfile=outdir + '/TaxaGeneIntegrity.rmsp.png', dpi=300, nfsize=30, tickstep=1)
    else:
        # Process alignments without spurious sequence removal
        dq = lambda x: ReadSeq(x, exclude=excl, outdir=outdir, mingenes=1)
        AlignSeries = list(map(dq, SelectedLoci['fastapath']))
        AlignSeries = pd.Series(AlignSeries, index=[x.name for x in AlignSeries])
        AlignSeries = AlignSeries[AlignSeries.apply(len) >= 10]  # Require at least 10 taxa
        
        # Concatenate and analyze
        AlignLength = ConcatFas2Phy(AlignSeries, outfile=outdir + '/ConcatenatedMatrix', 
                                   partitionfile=True, returnaspyobj=True)
        AlignLength.to_excel(outdir + '/MultipleSequenceAlignmentsLogLength.xlsx')
        HeatMap(AlignLength, outfile=outdir + '/TaxaGeneIntegrity.png', dpi=300, nfsize=30, tickstep=1)
    

