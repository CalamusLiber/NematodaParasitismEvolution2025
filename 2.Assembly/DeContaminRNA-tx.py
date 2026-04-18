#!/To/Your/Directory/anaconda3/bin/python3.9

"""
Decontamination script for removing vector and taxonomic contamination from assembled RNA sequences.

This script supports two main modes:
- Vector removal only when --valid is unset
- Taxonomic filtering when --valid is provided

It reads an input FASTA assembly and a BLAST/DIAMOND contamination table,
then writes filtered sequences to an output FASTA file.
"""

# Load required libraries.
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import pandas as pd
import argparse, os
from math import ceil

# Set the email address required by NCBI Entrez.
Entrez.email = 'lianges.luex@gmail.com'

def qsmaker(i):
    """Extract taxonomic IDs from a DIAMOND/BLAST row and attach them to a contig ID."""
    qs = pd.Series(str(Dec.iloc[i]['staxids']).split(';')).apply(lambda x: x if x.isdigit() else None).dropna()
    qs.index = [Dec.index[i]] * len(qs)
    return qs


def idcorr(oldid):
    """Convert merged taxonomy IDs to their current valid ID using merged.dmp."""
    newid = TaxMerge.loc[oldid]['validid'] if oldid in TaxMerge.index else oldid
    return newid


def getHigherLineageID(taxid, rank=''):
    """Return the lineage IDs for a taxon up to the root, optionally filtering by rank."""
    LineageID = [taxid]
    while LineageID[-1] != 1:
        LineageID.append(TaxNodes.loc[LineageID[-1]]['parent'])
    LineageID = pd.Series(LineageID, index=TaxNodes.loc[LineageID]['rank'])
    if rank == '':
        return LineageID
    elif rank in LineageID.index:
        return LineageID.iloc[LineageID.index.tolist().index(rank):]
    else:
        return None

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
    newpathDict = {k: [i - slice[0] if i >= slice[0] else i for i in v] for k, v in oripathDict.items() if v[0] < slice[1] and v[1] >= slice[0]}
    newpathDict = {k: [slice[1] - slice[0] - 1 if i >= slice[1] - slice[0] else i for i in v] for k, v in newpathDict.items()}
    newpathDict = {k: '-'.join([str(i) for i in v]) for k, v in newpathDict.items()}
    newpath = ' path=[' + ' '.join([':'.join([k, v]) for k, v in newpathDict.items()]) + ']'
    return newpath

def RmVec(fas, dec, minlen=300, maxgc=95):
    """Remove vector-contaminated fragments from contigs and return filtered sequences."""
    VecMatch = set(dec.index)
    NewFas = []
    for contigi in fas:
        if contigi.id in VecMatch:
            VecRange = pd.DataFrame(dec.loc[contigi.id])
            if VecRange.shape == (2, 1):
                VecRange = VecRange.T
            else:
                VecRange = merge(VecRange)
            ContigLen = len(contigi.seq)
            DeVecInd = split(ContigLen, VecRange)
            DeVecInd = [i for i in DeVecInd if i[2] >= minlen]
            OriPath = contigi.description.split('path=')[1]
            NewSeq = [contigi.seq[i[0]:i[1]] for i in DeVecInd]
            NewDesc = ['len=' + str(DeVecInd[i][2]) + mappath(OriPath, DeVecInd[i]) for i in range(len(DeVecInd)) if GC(NewSeq[i]) <= maxgc]
            NewSeq = [i for i in NewSeq if GC(i) <= maxgc]
            # Create new fragment IDs for split contigs and preserve lengths/paths in the description.
            NewId = [contigi.id + '_x' + str(i + 1) if len(NewSeq) > 1 else contigi.id for i in range(len(NewSeq))]
            NewRec = [SeqRecord(NewSeq[i], id=NewId[i], description=NewDesc[i]) for i in range(len(NewSeq))]
            NewFas.extend(NewRec)
        elif len(contigi.seq) >= minlen and GC(contigi.seq) <= maxgc:
            # Keep contigs with no vector hits that satisfy length and GC filters.
            NewFas.append(contigi)
    return NewFas

def whitelist_remote(qsid, taxref, stepmax=100, std=1):
    """Validate taxonomic IDs using remote NCBI taxonomy queries."""
    SidUni = qsid.drop_duplicates()
    stp1 = [x * stepmax for x in range(ceil(len(SidUni) / stepmax))]
    stp2 = stp1[1:] + [len(SidUni)]
    hit = []
    for i in range(len(stp1)):
        handle = Entrez.efetch(db="taxonomy", id=SidUni.iloc[stp1[i]:stp2[i]].tolist())
        TaxInfo = Entrez.read(handle)
        handle.close()
        hit.extend([taxref in t['Lineage'].split('; ') for t in TaxInfo])
    qsid = qsid.isin(SidUni[hit])
    if std == 0:
        vt = lambda x: any(x)  # relaxed sifting
    elif std == 1:
        vt = lambda x: sum(x) / len(x) >= 0.5  # moderate sifting
    elif std == 2:
        vt = lambda x: all(x)  # strict sifting
    QSidvote = qsid.groupby([qsid.index]).apply(vt)
    whitelist = list(QSidvote[QSidvote].index)
    return whitelist

# The local validation function can reduce false positives by distinguishing homonymous taxa across kingdoms.
def whitelist_local(qsid, taxref, std=1):
    """Validate taxonomic IDs against a local NCBI taxonomy dump."""
    qtaxrefId = TaxNames[(TaxNames['nametxt'] == taxref) & (TaxNames['nameclass'] == 'scientific name')].index.values[0]
    qhigertaxId = getHigherLineageID(qtaxrefId)
    SidUni = qsid.drop_duplicates()
    hit = SidUni.apply(lambda x: all(qhigertaxId.isin(getHigherLineageID(x))))
    qsid = qsid.isin(SidUni[hit])
    if std == 0:
        vt = lambda x: any(x)  # relaxed sifting
    elif std == 1:
        vt = lambda x: sum(x) / len(x) >= 0.5  # moderate sifting
    elif std == 2:
        vt = lambda x: all(x)  # strict sifting
    QSidvote = qsid.groupby([qsid.index]).apply(vt)
    whitelist = list(QSidvote[QSidvote].index)
    return whitelist

def retagcontigs(fas):
    """Restore original contig names when only one split fragment remains."""
    Xtag = [c.id for c in fas if c.id.find('_x') != -1]
    Otag = [d.split('_x')[0] for d in Xtag]
    SiXtag = {t for t in Xtag if Otag.count(t.split('_x')[0]) == 1}
    for f in range(len(fas)):
        if fas[f].id in SiXtag:
            fas[f].id = fas[f].id.split('_x')[0]
    return fas


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Decontaminate assembled RNA sequences using vector and taxonomic filters.")
    parser.add_argument('--fasta', '-f',
                        type=str,
                        help='Input FASTA file containing sequences to decontaminate.')
    parser.add_argument('--csv', '-c',
                        type=str,
                        help='Tab-delimited CSV file with contaminated fragment coordinates from BLAST.')
    parser.add_argument('--out', '-o',
                        type=str,
                        help='Output directory or filename for the decontaminated sequences.',
                        default=".")
    parser.add_argument('--keep', '-k',
                        action='store_true',
                        help='Keep internally trimmed fragments in the output along with intact contigs.')
    parser.add_argument('--splitout', '-s',
                        action='store_true',
                        help='Also write trimmed fragments to a separate Split.fasta file.')
    parser.add_argument('--sensitive', '-t',
                        type=str,
                        help='Whitelist sensitivity: 0=relaxed, 1=moderate, 2=strict.',
                        default="1")
    parser.add_argument('--evaluecutoff', '-e',
                        type=str,
                        help='Ignore BLAST hits with E-value larger than this cutoff (default 1e-20).',
                        default="1e-20")
    parser.add_argument('--valid', '-v',
                        type=str,
                        help='NCBI taxon name for validation; if set, perform taxonomic filtering instead of only vector removal.',
                        default="")
    parser.add_argument('--db', '-d',
                        type=str,
                        help='Taxonomy source: "remote" for online NCBI queries, or local taxdump directory path.',
                        default="remote")
    args = parser.parse_args()

    fastafile = args.fasta.replace('\\', '/')
    decontfile = args.csv
    if args.out == '':
        outfile = fastafile.replace('.fasta', ".Sifted.fasta")
    else:
        outfile = args.out + "/" + fastafile.split("/")[-1].replace('.fasta', ".Sifted.fasta")

    # Load BLAST contamination coordinates and filter by E-value.
    Dec = pd.read_csv(decontfile, index_col=0, names=['contigid', 'sseqid', 'qstart', 'qend', 'score', 'evalue', 'pident', 'ppos', 'staxids'], sep="\t")
    Dec = Dec[Dec['evalue'] <= float(args.evaluecutoff)]

    # Load input sequences.
    Fas = list(SeqIO.parse(fastafile, "fasta"))

    if args.valid == "":
        # Vector-removal mode: remove vector-contaminated regions only.
        Dec = Dec[['qstart', 'qend']]
        Dec.sort_index(inplace=True)
        FasDeVec = RmVec(Fas, Dec, minlen=300, maxgc=95)
        if args.keep:
            SeqIO.write(FasDeVec, outfile, "fasta")
            if args.splitout:
                FasSplit = [i for i in FasDeVec if i.id.find('_x') != -1]
                SeqIO.write(FasSplit, outfile.replace('.Sifted.fasta', ".Split.fasta"), "fasta")
        elif args.splitout:
            FasSplit = [i for i in FasDeVec if i.id.find('_x') != -1]
            FasDeVec = [i for i in FasDeVec if i.id.find('_x') == -1]
            SeqIO.write(FasDeVec, outfile, "fasta")
            SeqIO.write(FasSplit, outfile.replace('.Sifted.fasta', ".Split.fasta"), "fasta")
        else:
            FasDeVec = [i for i in FasDeVec if i.id.find('_x') == -1]
            SeqIO.write(FasDeVec, outfile, "fasta")
    else:
        # Taxonomic filtering mode: keep only contigs assigned to the designated taxon.
        qs = list(map(qsmaker, range(len(Dec.index))))
        QSid = pd.concat(qs)
        QSid.dropna(inplace=True)
        qreftax = args.valid.capitalize()
        if args.db.lower() == "remote":
            # Query NCBI Taxonomy remotely.
            QSid = QSid.astype(str)
            whl = whitelist_remote(qsid=QSid, taxref=qreftax, std=int(args.sensitive), stepmax=50)
        elif all([os.path.exists(args.db + '/nodes.dmp'), os.path.exists(args.db + '/names.dmp'), os.path.exists(args.db + '/merged.dmp'), os.path.exists(args.db + '/delnodes.dmp')]):
            # Use a local taxonomy dump if available.
            TaxNodes = pd.read_csv(args.db + '/nodes.dmp', names=['taxid', 'parent', 'rank'], index_col=0, header=None, usecols=[0,2,4], sep="\t")
            TaxNames = pd.read_csv(args.db + '/names.dmp', names=['taxid', 'nametxt', 'uniquename', 'nameclass'], index_col=0, header=None, usecols=[0,2,4,6], sep="\t")
            TaxMerge = pd.read_csv(args.db + '/merged.dmp', names=['oldid', 'validid'], index_col=0, header=None, usecols=[0,2], sep="\t")
            TaxDel = pd.read_csv(args.db + '/delnodes.dmp', index_col=None, header=None, usecols=[0], sep="\t").iloc[:,0]
            QSid = QSid.astype(int)
            QSid = QSid.apply(idcorr)
            QSid = QSid[~QSid.isin(TaxDel)]
            whl = whitelist_local(qsid=QSid, taxref=qreftax, std=int(args.sensitive))
        else:
            print('Wrong input for the argument of "--db".\n')
            exit()
        FasDeCon = [x for x in Fas if x.id in whl]
        FasDeCon = retagcontigs(FasDeCon)
        SeqIO.write(FasDeCon, outfile.replace('.Sifted.fasta', ".DeContX.fasta"), "fasta")
