#!/To/Your/Directory/anaconda3/bin/python3.9

"""
UniteBUSCOs.py - Consolidate BUSCO Results Across Multiple Taxa

This script processes BUSCO (Benchmarking Universal Single-Copy Orthologs) results
from multiple taxa/species. It consolidates protein sequences for each BUSCO gene
by selecting the longest isoform across all available predictions for each taxon.

The script performs the following operations:
1. Maps directories to taxonomic groups based on naming patterns
2. Identifies BUSCO gene FASTA files in each directory
3. For each BUSCO gene, selects the longest protein sequence per taxon
4. Generates unified FASTA files for each BUSCO gene
5. Creates statistical summaries and visualizations of gene presence/absence

Input: Directory containing subdirectories with BUSCO protein predictions (*.faa files)
Output: Unified gene FASTA files, Excel summaries, and heatmap visualization

Dependencies: Biopython, matplotlib, pandas, numpy
"""

import os, re
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib import gridspec as gs
import numpy as np
import pandas as pd
import argparse

def TaxaDirMap(dir):
    """
    Map directory names to taxonomic groups.
    
    Parses directory names to extract taxon identifiers by splitting on common
    accession prefixes (SRS, ERP, ERS, GCA, GCF, etc.) and groups directories
    belonging to the same taxon.
    
    Args:
        dir (pd.Series): Series of directory names
        
    Returns:
        pd.Series: Series with taxon names as index and lists of directories as values
    """
    # Extract taxon names by splitting on accession prefixes
    taxa = dir.str.split(r'(_SRS|_ERP|_ERS|_GCA|_GCF|_UP|.MG|_PRJEB|_NC|_NEB|_N1|_N2)', 
                        expand=True, n=1).iloc[:,0].drop_duplicates()
    taxa.reset_index(drop=True, inplace=True)    
    
    # Group directories by taxon
    tdlist = taxa.apply(lambda x: dir[dir.str.find(x)==0])
    tdlist.index = taxa.values
    return tdlist

def DirFastaMap(dir, prefix='.'):
    """
    Map directories to their contained FASTA files.
    
    For each directory, identifies all .faa (protein FASTA) files.
    
    Args:
        dir (pd.Series): Series of directory names
        prefix (str): Base path prefix for directories
        
    Returns:
        pd.Series: Series with directory names as index and DataFrames of [dir, file] pairs as values
    """
    dflist = dir.apply(lambda x: pd.DataFrame([[x, i] for i in os.listdir(prefix+'/'+x) if i[-4:]=='.faa']))
    dflist.index = dir.values
    return dflist

def TaxaFastaTable(dirfasta, prefix='.'):
    """
    Create a comprehensive table of taxa, genes, and file paths.
    
    Consolidates all directory-FASTA mappings into a single table with
    taxon names, gene identifiers, and full file paths.
    
    Args:
        dirfasta (pd.Series): Series from DirFastaMap output
        prefix (str): Base path prefix
        
    Returns:
        pd.DataFrame: Table with columns 'Taxon', 'Gene', 'Path'
    """
    df = pd.concat(dirfasta.tolist(), ignore_index=True)
    
    # Extract taxon names from directory names
    taxa = df.iloc[:,0].str.split(r'(_SRS|_ERP|_ERS|_GCA|_GCF|_UP|.MG|_PRJEB|_NC|_NEB|_N1|_N2)', 
                                 expand=True, n=1).iloc[:,0]
    gene = df.iloc[:,1]
    
    # Construct full file paths
    fastapath = df.apply(lambda x: prefix + '/' + x[0] + '/' + x[1], axis=1)
    
    return pd.DataFrame({'Taxon':taxa, 'Gene':gene, 'Path':fastapath}, index=range(len(taxa)))

def selectLongest(file, id):
    """
    Select the longest protein sequence from multiple FASTA files.
    
    Given a list of FASTA files for the same BUSCO gene from different
    predictions, selects the one with the longest sequence.
    
    Args:
        file (list): List of file paths to FASTA files
        id (str): Identifier to assign to the selected sequence
        
    Returns:
        SeqRecord: Biopython SeqRecord with the longest sequence
    """
    seqlen = 0
    for i in file:
        fa = SeqIO.read(i, "fasta")
        if len(fa) > seqlen:
            fasta = fa
            fasta.id = id  # Set sequence ID to taxon name
            fasta.name = ''
            fasta.description = ''
            seqlen = len(fasta)
    return fasta

def writeSeq(gene):
    """
    Process a single BUSCO gene across all taxa.
    
    For a given gene, collects all available sequences from different taxa,
    selects the longest one for each taxon, and writes a unified FASTA file.
    
    Args:
        gene (str): BUSCO gene identifier
        
    Returns:
        pd.DataFrame: Statistics table for this gene across taxa
    """
    # Get taxa that have this gene
    fp = TaxaGeneUni.loc[TaxaGeneUni.loc[:,'Gene']==gene,:]
    
    # Select longest sequence for each taxon
    genefasta = fp.apply(lambda x: selectLongest(
        file=TaxaFasta.loc[(TaxaFasta.loc[:,'Taxon']==x['Taxon']) & 
                           (TaxaFasta.loc[:,'Gene']==gene),'Path'].tolist(), 
        id=x['Taxon']), axis=1)
    
    # Write unified FASTA file
    SeqIO.write(genefasta, OutDir+'/BSC_'+gene.replace('at6231',''), "fasta")
    
    # Create statistics table
    res = fp.copy()
    res.loc[:,'Max.Length'] = genefasta.apply(len)
    return res

def TaxaGeneMatrix(genestattable):
    """
    Create a matrix of taxa vs genes with sequence lengths.
    
    Transforms the gene statistics table into a matrix format where
    rows are taxa and columns are genes, with cell values being
    sequence lengths (NaN for missing genes).
    
    Args:
        genestattable (pd.DataFrame): Gene statistics table from writeSeq
        
    Returns:
        pd.DataFrame: Taxa-gene matrix with sequence lengths
    """
    genestattable.index = genestattable.loc[:,'Taxon']
    genuni = genestattable.loc[:,'Gene'].drop_duplicates()
    
    # Create matrix by pivoting gene lengths
    tgmat = list(map(lambda x: genestattable.loc[genestattable.loc[:,'Gene']==x, 'Max.Length'], genuni))
    tgmat = pd.concat(tgmat, sort=True, axis=1)
    tgmat.columns = genuni
    tgmat.sort_index(axis=1, inplace=True)
    
    # Add total length column
    totlen = tgmat.apply(lambda x: sum(x.dropna()), axis=1)
    totlen.name = 'Total.Length'
    tgmat = pd.concat([totlen, tgmat], sort=True, axis=1)
    tgmat.reset_index(inplace=True)    
    return tgmat

def HeatMap(taxastat, outdir, dpi=300, nfsize=50, tickstep=10):
    """
    Generate a three-panel heatmap visualization of BUSCO results.
    
    Creates a figure with three subplots:
    1. Log-transformed sequence lengths (heatmap)
    2. Gene presence/absence (binary heatmap)  
    3. Bar chart of missing genes per BUSCO
    
    Args:
        taxastat (pd.DataFrame): Taxa-gene statistics matrix
        outdir (str): Output directory for the plot
        dpi (int): Resolution for saved image
        nfsize (int): Figure size multiplier
        tickstep (int): Step size for axis tick labels
    """
    # Prepare data for plotting
    ts_log = taxastat.iloc[:-2,2:].apply(np.log).fillna(0)  # Log lengths, fill NA with 0
    ts_bin = taxastat.iloc[:-2,2:].isnull()  # Binary presence/absence
    
    # Set up tick positions
    xlabpos = np.arange(len(ts_log.columns), step=tickstep*3)
    ylabpos = np.arange(len(ts_log.index), step=tickstep)
    
    # Create figure with three subplots
    fig = plt.figure(figsize=(3*nfsize,2*nfsize))
    gsp = gs.GridSpec(3, 1, height_ratios=[1, 1, 1])
    
    # Top plot: Log sequence lengths
    axes0 = plt.subplot(gsp[0])
    htp0 = axes0.imshow(ts_log, cmap=plt.cm.jet, aspect='auto')
    axes0.set_title('Log-length of gene vs. taxon', fontdict={'fontsize':25})
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
    axes2.bar(taxastat.iloc[-1,2:].index, taxastat.iloc[-1,2:], width=1, color='orange')
    axes2.margins(x=0)
    axes2.set_xlabel('BUSCO gene', fontdict={'fontsize':18})
    axes2.set_ylabel('Missing', fontdict={'fontsize':18})
    axes2.set_title('Count of missing genes', fontdict={'fontsize':25})
    axes2.set_xticks(xlabpos)
    axes2.set_xticklabels(labels=ts_bin.columns[xlabpos].tolist(), fontdict={'fontsize':3})
    plt.setp(axes2.get_xticklabels(), rotation=270)
    
    # Finalize and save plot
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    fig.savefig(outdir+'/'+'TaxaGeneCount.png', dpi=dpi)

if __name__ == "__main__":
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Consolidate BUSCO results across multiple taxa by selecting longest sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python UniteBUSCOs.py -i busco_results/ -o united_buscos/
  python UniteBUSCOs.py -i ./busco_output -o ./consolidated -s
        """
    )
    parser.add_argument('--indir', '-i',
                        type=str,
                        required=True,
                        help='Input directory containing BUSCO result subdirectories')
    parser.add_argument('--outdir', '-o',
                        type=str,
                        help='Output directory for consolidated results (default: current directory)',
                        default=".")
    parser.add_argument('--sortmissing', '-s',
                        action='store_true',
                        help='Sort output by number of missing genes (ascending) instead of alphabetical')
    
    args = parser.parse_args()

    # Process input/output paths
    InDir = args.indir.replace('\\','/')
    OutDir = args.outdir.replace('\\','/')
    CurDir = os.getcwd().replace('\\','/')
    
    # Convert to absolute paths if relative
    if not os.path.isabs(InDir):
        InDir = CurDir + '/' + InDir.replace('./','')
    if not os.path.isabs(OutDir):
        OutDir = CurDir + '/' + OutDir.replace('./','')
    
    # Create output directory if it doesn't exist
    if not os.path.exists(OutDir):
        os.makedirs(OutDir, mode=0o777, exist_ok=True)    
    
    # Initialize data collection
    DirList = os.listdir(InDir)
    DirList = pd.Series([i for i in DirList if os.path.isdir(InDir+'/'+i)])
    
    # Map directories to taxa and identify FASTA files
    TaxaDir = TaxaDirMap(DirList)
    DirFasta = DirFastaMap(DirList, prefix=InDir)
    TaxaFasta = TaxaFastaTable(DirFasta, prefix=InDir)
    
    # Count genes per taxon and get unique combinations
    GeneCount = TaxaFasta.groupby(['Taxon','Gene']).count().reset_index()
    TaxaGeneUni = TaxaFasta[['Taxon','Gene']].drop_duplicates(inplace=False, ignore_index=True)
    GeneUni = TaxaGeneUni['Gene'].drop_duplicates(inplace=False)

    # Generate unified gene FASTA files and collect statistics
    print("Processing BUSCO genes...")
    GeneLenght = list(map(writeSeq, GeneUni))
    GeneLenght = pd.concat(GeneLenght, ignore_index=True)
    GeneStat = pd.merge(GeneCount, GeneLenght, on=['Taxon', 'Gene'], sort=True)
    GeneStat['Gene'] = GeneStat['Gene'].apply(lambda x: 'BSC_'+x.replace('at6231.faa',''))
    GeneStat.columns = ['Taxon', 'Gene', 'Source', 'Max.Length']
    del(GeneLenght)

    # Generate taxa-gene summary statistics
    TaxaStat1 = TaxaGeneUni.groupby(['Taxon']).count().reset_index()
    TaxaStat1.sort_values(by='Taxon', ascending=True, inplace=True)
    TaxaStat1.columns = ['Taxon', 'Total.Genes']
    TaxaStat2 = TaxaGeneMatrix(GeneStat)
    TaxaStat = pd.merge(TaxaStat1, TaxaStat2, on=['Taxon'], sort=True)
    TaxaStat.index = TaxaStat.loc[:, 'Taxon']
    TaxaStat = TaxaStat.iloc[:, 1:]
    
    # Add summary statistics rows
    TSavg = pd.Series(TaxaStat.median(axis=0), name='MEDIAN').astype(int)
    TSmissing = pd.Series(TaxaStat.isnull().sum(axis=0), name='MISSING').astype(int)
    TSmissing[0] = sum(TSmissing[2:])  # Total missing across all genes
    TaxaStat.loc['AVERAGE',:] = TSavg
    TaxaStat.loc['MISSING',:] = TSmissing
    
    # Optional sorting by missing genes
    if args.sortmissing:
        colindex = TaxaStat.columns[:2].tolist() + TaxaStat.iloc[-1, 2:].sort_values(ascending=True).index.tolist()
        rowindex = TaxaStat.iloc[:-2, 0].sort_values(ascending=False).index.tolist() + TaxaStat.index[-2:].tolist()
        TaxaStat = TaxaStat.loc[rowindex, colindex]
    
    # Clean up temporary variables
    del(TaxaStat1, TaxaStat2, TSavg, TSmissing)

    # Write output files
    GeneStat.to_excel(OutDir+'/'+'UnitedGeneStat.xlsx')
    TaxaStat.to_excel(OutDir+'/'+'TaxaGeneCount.xlsx')
    
    # Print summary statistics
    print('In summary, the UniteBUSCOs processing has treated {DL} protein prediction files and extracted {NG} genes for {NT} species. In the output {NT} * {NG} data matrix, there are {NM} missing loci ({PM}%).\n'.format(
        DL=str(len(DirList)), 
        NG=str(len(GeneUni)), 
        NT=str(len(TaxaStat.index)-2), 
        NM=str(TaxaStat.isna().sum().sum()), 
        PM=str(round(TaxaStat.isna().sum().sum()*100/(len(TaxaStat.index)-2)/len(GeneUni),2))))
    
    # Generate visualization
    print("Creating heatmap visualization...")
    HeatMap(taxastat=TaxaStat, outdir=OutDir, nfsize=30, tickstep=1)
