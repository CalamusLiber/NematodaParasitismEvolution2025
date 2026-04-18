#!/To/Your/Directory/anaconda3/bin/python3.9

"""
GffRepeatMasker.py - Sequence Masking Tool for Repeat Elements

This script masks repetitive or contaminated regions in FASTA sequences based on
GFF/RepeatMasker output coordinates. It supports both hard masking (replacing with N's)
and soft masking (converting to lowercase).

The script reads a FASTA file containing sequences and a GFF/CSV file with repeat
coordinates, then applies the specified masking strategy to regions identified
as repetitive or contaminated.

Usage:
    python GffRepeatMasker.py --fasta input.fasta --gff repeats.gff --masktype hard

Dependencies: Biopython, pandas
"""

# Load required libraries
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import pandas as pd
import argparse

def toN(seq, gffline):
    """
    Replace a sequence region with N characters (hard masking).
    
    Args:
        seq (MutableSeq): The mutable sequence to modify
        gffline (pd.Series): Row from GFF dataframe with 'qstart' and 'qend' columns
    """
    s = gffline['qstart']-1  # Convert to 0-based indexing
    e = gffline['qend']    
    seq[s:e] = 'N' * (e - s)  # Replace region with N's

def toLower(seq, gffline):
    """
    Convert a sequence region to lowercase (soft masking).
    
    Args:
        seq (MutableSeq): The mutable sequence to modify
        gffline (pd.Series): Row from GFF dataframe with 'qstart' and 'qend' columns
    """
    s = gffline['qstart']-1  # Convert to 0-based indexing
    e = gffline['qend']
    seq[s:e] = seq[s:e].lower()  # Convert region to lowercase

def HardMask(contig):
    """
    Apply hard masking to a single contig/sequence.
    
    Replaces identified repeat regions with N characters. If multiple repeat
    regions exist for this contig, all are masked.
    
    Args:
        contig (SeqRecord): Biopython SeqRecord object to be masked
    """
    seqid = contig.id
    contig.seq = contig.seq.upper()  # Ensure sequence is uppercase first
    
    if seqid in Gff.index:
        seqnew = MutableSeq(contig.seq)
        repregion = Gff.loc[seqid]
        
        # Handle single vs multiple repeat regions
        if repregion.shape==(2,):
            toN(seqnew, repregion)  # Single region
        else:
            repregion.apply(lambda x: toN(seqnew, x), axis=1)  # Multiple regions
        
        contig.seq = Seq(seqnew)
    else:
        pass  # No repeats found for this contig

def SoftMask(contig):
    """
    Apply soft masking to a single contig/sequence.
    
    Converts identified repeat regions to lowercase. If multiple repeat
    regions exist for this contig, all are masked.
    
    Args:
        contig (SeqRecord): Biopython SeqRecord object to be masked
    """
    seqid = contig.id
    contig.seq = contig.seq.upper()  # Ensure sequence is uppercase first
    
    if seqid in Gff.index:
        seqnew = MutableSeq(contig.seq)
        repregion = Gff.loc[seqid]
        
        # Handle single vs multiple repeat regions
        if repregion.shape==(2,):
            toLower(seqnew, repregion)  # Single region
        else:
            repregion.apply(lambda x: toLower(seqnew, x), axis=1)  # Multiple regions
        
        contig.seq = Seq(seqnew)
    else:
        pass  # No repeats found for this contig

if __name__ == "__main__":
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Mask repetitive or contaminated regions in FASTA sequences based on GFF coordinates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python GffRepeatMasker.py -f genome.fasta -g repeats.gff -t hard -o masked/
  python GffRepeatMasker.py -f transcripts.fasta -g contamination.csv -t both -s 1
        """
    )
    parser.add_argument('--fasta', '-f',
                        type=str,
                        required=True,
                        help='Input FASTA file containing sequences to be masked')
    parser.add_argument('--gff', '-c',
                        type=str,
                        required=True,
                        help='GFF/CSV file with repeat/contamination coordinates (tab-separated)')
    parser.add_argument('--out', '-o',
                        type=str,
                        help='Output directory for masked sequences (default: current directory)',
                        default=".")
    parser.add_argument('--skip', '-s',
                        type=str,
                        help='Number of header lines to skip when reading GFF file (default: 0)',
                        default="0")
    parser.add_argument('--masktype', '-t',
                        type=str,
                        choices=['hard', 'soft', 'both'],
                        help='Masking type: hard (replace with N), soft (lowercase), both (create both outputs)',
                        default="soft")
    
    args = parser.parse_args()

    # Process file paths (handle Windows backslashes)
    fastafile = args.fasta.replace('\\','/')
    fmt = fastafile.split('.')[-1]  # Extract file extension
    gfffile = args.gff.replace('\\','/')
    
    # Determine output directory
    if args.out=='':
        outprefix = fastafile.replace(fastafile.split("/")[-1], "")
    else:
        outprefix = args.out.replace('\\','/') + '/'
    
    # Load GFF data: sequence ID, start position, end position
    # Assumes tab-separated format with columns: seqid, ..., qstart, qend
    Gff = pd.read_csv(gfffile, names=['sseqid', 'qstart', 'qend'], 
                     index_col=0, header=None, usecols=[0,3,4], 
                     skiprows=int(args.skip), sep="\t")
    
    # Load FASTA sequences into pandas Series for easy processing
    Fas = pd.Series(SeqIO.parse(fastafile, "fasta"))
    
    # Apply masking based on user selection
    if args.masktype.lower()=='hard':
        Fas.apply(HardMask)
        SeqIO.write(Fas, outprefix + fastafile.split("/")[-1].replace('.'+fmt, ".hardmask."+fmt), "fasta")
    elif args.masktype.lower()=='soft':
        Fas.apply(SoftMask)
        SeqIO.write(Fas, outprefix + fastafile.split("/")[-1].replace('.'+fmt, ".softmask."+fmt), "fasta")
    elif args.masktype.lower()=='both':
        # Create both soft and hard masked versions
        Fas.apply(SoftMask)
        SeqIO.write(Fas, outprefix + fastafile.split("/")[-1].replace('.'+fmt, ".softmask."+fmt), "fasta")
        Fas.apply(HardMask)  # Apply hard masking to the already soft-masked sequences
        SeqIO.write(Fas, outprefix + fastafile.split("/")[-1].replace('.'+fmt, ".hardmask."+fmt), "fasta")
    else:
        print("Error: Invalid --masktype setting. Please use 'hard', 'soft', or 'both'.")
        quit()
