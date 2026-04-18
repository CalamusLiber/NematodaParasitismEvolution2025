#!/To/Your/Directory/miniconda3/bin/python3.9

"""
Multi-threaded vector decontamination script for assembled RNA sequences.

This script removes BLAST-identified vector-contaminated regions from assembled contigs
using parallel processing. It reads a FASTA file and a BLAST alignment table,
then outputs clean sequences with vector regions trimmed.
"""

# Load required libraries.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
from multiprocessing.dummy import Pool as ThreadPool

# Example usage (commented out):
# seq_tri = list(SeqIO.parse("Euperipatoides_rowelli_SRS6182028_trinity.Trinity.fasta", "fasta"))
# DeCont = pd.read_csv("Euperipatoides_rowelli_SRS6182028_Trinity.DeCont", index_col=0, names=['sub', 'identity', 'alignlen', 'mismatched', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

def merge(Rg):
    """Merge overlapping or adjacent vector hit intervals into non-redundant ranges."""
    Rg.sort_values(by=['qstart', 'qend'], axis=0, inplace=True)
    for i in range(1, len(Rg)):
        if Rg.iloc[i, 0] <= Rg.iloc[i - 1, 1]:
            Rg.iloc[[i - 1, i], 1] = max(Rg.iloc[[i - 1, i], 1])
            Rg.iloc[i, 0] = Rg.iloc[i - 1, 0]
    Rg.drop_duplicates(inplace=True)
    Rg.sort_values(by=['qstart', 'qend'], inplace=True)
    return Rg


def split(seqlen, rg):
    """Return contiguous sequence segments after removing the matched regions."""
    ht = pd.DataFrame([[-1, 0], [seqlen + 1, -1]], columns=['qstart', 'qend'])
    rg = pd.concat([ht, rg])
    rg.sort_values(by=['qstart', 'qend'], axis=0, inplace=True)
    kp = []
    for i in range(len(rg) - 1):
        d = rg.iloc[i + 1, 0] - rg.iloc[i, 1] - 1
        if d > 0:
            kp.append([int(rg.iloc[i, 1]), int(rg.iloc[i + 1, 0] - 1), int(d)])
    return kp

def mappath(oripath, slice):
    """Translate original path annotations to coordinates for a sliced contig fragment."""
    oripathDict = {i.split(':')[0]: list(map(int, i.split(':')[1].split('-'))) for i in oripath[1:-1].split()}
    newpathDict = {k: [i - slice[0] if i >= slice[0] else i for i in v] for k, v in oripathDict.items() if slice[0] <= v[0] < slice[1] or slice[0] <= v[1] < slice[1]}
    newpathDict = {k: [slice[1] - slice[0] - 1 if i >= slice[1] - slice[0] else i for i in v] for k, v in newpathDict.items()}
    newpathDict = {k: '-'.join([str(i) for i in v]) for k, v in newpathDict.items()}
    newpath = ' path=[' + ' '.join([':'.join([k, v]) for k, v in newpathDict.items()]) + ']'
    return newpath

def RmVec(contigi, dec, vecm, minlen=300):
    """Remove vector-contaminated regions from a single contig and return cleaned fragments."""
    NewRec = []
    if contigi.id in vecm:
        VecRange = pd.DataFrame(dec.loc[contigi.id])
        if VecRange.shape == (2, 1):
            VecRange = VecRange.T
        else:
            VecRange = merge(VecRange)
        ContigLen = len(contigi.seq)
        DeVecInd = split(ContigLen, VecRange)
        DeVecInd = [i for i in DeVecInd if i[2] >= minlen]
        OriPath = contigi.description.split('path=')[1]
        
        # Create new fragments for each valid cleaned segment.
        for i in range(len(DeVecInd)):
            NewSeq = contigi.seq[DeVecInd[i][0]:DeVecInd[i][1]]
            NewPath = mappath(OriPath, DeVecInd[i])
            NewId = contigi.id + '_x' + str(i + 1)
            NewDesc = NewId + ' len=' + str(DeVecInd[i][2]) + NewPath
            NewRec.append(SeqRecord(NewSeq, id=NewId, name=NewId, description=NewDesc))
    elif len(contigi.seq) >= minlen:
        # Keep contigs that have no vector hits and meet the minimum length.
        NewRec = contigi
    return NewRec

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove BLAST-identified vector regions from assembled RNA contigs using multi-threaded processing.")
    parser.add_argument('--fasta', '-f',
                        type=str,
                        help='Input FASTA file containing sequences to decontaminate.')
    parser.add_argument('--csv', '-c',
                        type=str,
                        help='Tab-delimited CSV file with vector contamination coordinates from BLAST (columns: query_id, qstart, qend).')
    parser.add_argument('--out', '-o',
                        type=str,
                        help='Output directory for the decontaminated sequences.',
                        default=".")
    args = parser.parse_args()

    fastafile = args.fasta
    decontfile = args.csv
    outfile = args.out + "/" + fastafile.split("/")[-1].replace('.fasta', ".Sifted.fasta")
    decontheader = ['', 'qstart', 'qend']

    # Load input sequences and contamination coordinates.
    Fas = list(SeqIO.parse(fastafile, "fasta"))
    Dec = pd.read_csv(decontfile, index_col=0, names=decontheader, usecols=[0, 6, 7])
    Dec.sort_index(inplace=True)
    VecMatch = set(Dec.index)

    # Apply vector removal in parallel using 8 threads.
    pool = ThreadPool(8)
    FasDeVecX = [pool.apply(RmVec, args=(contigi, Dec, VecMatch, 300)) for contigi in Fas]
    pool.close()
    
    # Flatten the list of results and write to output file.
    FasDeVec = []
    for f in FasDeVecX:
        FasDeVec.extend(f)

    SeqIO.write(FasDeVec, outfile, "fasta")
