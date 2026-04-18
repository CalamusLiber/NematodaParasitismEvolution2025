#!/To/Your/Directory/anaconda3/bin/python3.9

"""
Generate pie charts visualizing taxonomic composition of decontaminated transcripts.

This script analyzes NCBI taxonomy assignments for sequence hits and creates pie charts
showing the distribution of transcript origins across major taxonomic groups.
Supports both remote NCBI Taxonomy queries and local taxonomy database lookups.
"""

# Load required libraries.
import matplotlib.pyplot as plt
import pandas as pd
import argparse, os
from math import ceil
from Bio import Entrez

# Set the email address required by NCBI Entrez API.
Entrez.email = 'lianges.luex@gmail.com'

def qsmaker(i):
    """Extract taxonomic IDs from a DIAMOND/BLAST row and associate them with the query contig."""
    qs = pd.Series(str(Dec.iloc[i]['staxids']).split(';')).apply(lambda x: x if x.isdigit() else None).dropna()
    qs.index = [Dec.index[i]] * len(qs)
    return qs

def idcorr(oldid):
    """Convert merged taxonomy IDs to their current valid ID using merged.dmp."""
    newid = TaxMerge.loc[oldid]['validid'] if oldid in TaxMerge.index else oldid
    return newid


def getHigherLineage(taxid, rank=''):
    """Return the scientific names in a taxon's lineage, optionally filtered by rank."""
    LineageID = [taxid]
    while LineageID[-1] != 1:
        LineageID.append(TaxNodes.loc[LineageID[-1]]['parent'])
    LineageID = pd.Series(LineageID, index=TaxNodes.loc[LineageID]['rank'])
    if rank == '':
        return TaxNames.loc[LineageID]['nametxt']
    elif rank in LineageID.index:
        return TaxNames.loc[LineageID.iloc[LineageID.index.tolist().index(rank):]]['nametxt']
    else:
        return None

def qclass(GenName, rank, remote=False):
    """Retrieve the taxonomic name at a specified rank for a query organism.
    
    Can use either remote NCBI Taxonomy API or local taxonomy database.
    """
    if remote:
        qtaxrefId = Entrez.read(Entrez.esearch(db='taxonomy', term=GenName))['IdList'][0]
        qtaxrefLin = Entrez.read(Entrez.efetch(db='taxonomy', id=qtaxrefId))[0]['LineageEx']
        qreftax = [rk['ScientificName'] for rk in qtaxrefLin if rk['Rank'].lower() == rank][0]
    else:
        qtaxrefId = TaxNames[TaxNames['nametxt'] == GenName].index.values[0]
        qreftax = getHigherLineage(qtaxrefId, rank=rank).iloc[0]
    return qreftax

def freq(dec, target, stepmax=100, remote=False):
    """Analyze taxonomic composition and generate frequency distribution across taxa.
    
    Extracts taxonomy IDs from BLAST hits, retrieves their lineages (remote or local),
    classifies them into major taxonomic groups, and performs voting-based consensus
    assignment for each query sequence.
    """
    # Extract all taxonomic IDs from query sequences.
    qs = list(map(qsmaker, range(len(dec.index))))
    qsid = pd.concat(qs)
    qsid.dropna(inplace=True)
    
    if remote:
        # Query NCBI Taxonomy remotely in batches.
        SidUni = qsid.drop_duplicates().astype(str)
        stp1 = [x * stepmax for x in range(ceil(len(SidUni) / stepmax))]
        stp2 = stp1[1:] + [len(SidUni)]
        rk = []
        for i in range(len(stp1)):
            handle = Entrez.efetch(db="taxonomy", id=SidUni[stp1[i]:stp2[i]].tolist())
            TaxInfo = Entrez.read(handle)
            handle.close()
            rk.extend([[t['AkaTaxIds'][0] if 'AkaTaxIds' in t.keys() else t['TaxId'], t['Lineage']] for t in TaxInfo])
        rk = pd.DataFrame(rk).drop_duplicates()  # 2-column dataframe with no index
        rk.index = rk.iloc[:, 0]
        rk = pd.Series(rk.iloc[:, 1].apply(lambda x: x.split('; ')))  # series with index
    else:
        # Use local taxonomy database for faster queries.
        qsid = qsid.astype(int)
        qsid = qsid.apply(idcorr)
        qsid = qsid[~qsid.isin(TaxDel)]
        SidUni = qsid.drop_duplicates()
        rk = SidUni.apply(lambda x: getHigherLineage(x).tolist())  # series with index
        rk.index = SidUni.values
        rk.dropna(inplace=True)
        rk.drop_duplicates(inplace=True)
    
    # Filter query IDs to those with valid taxonomy information.
    qsid = qsid[qsid.isin(rk.index)].astype(int)
    
    # Define taxonomic groups and their shorthand labels for visualization.
    TAGclade = pd.Series([target, 'Arthropoda', 'Ecdysozoa', 'Spiralia', 'Protostomia', 'Deuterostomia', 'Xenacoelomorpha', 'Bilateria', 'Cnidaria', 'Ctenophora', 'Placozoa', 'Porifera', 'Viridiplantae', 'Fungi', 'Bacteria', 'Archaea', 'Viruses'], index=[target, 'Arthr', 'Ecdys+', 'Spira', 'Proto+', 'Deute', 'Xenac', 'Bilat+', 'Cnida', 'Cteno', 'Placo', 'Porif', 'Virid', 'Fungi', 'Bacte', 'Archa', 'Virus'])
    
    def tagrk(x):
        """Map a subject taxonomy ID to the most specific major taxonomic group."""
        tagin = TAGclade.isin(rk[x])
        if any(tagin):
            return TAGclade.index[tagin][0]
        else:
            return 'others'
    
    # Assign each subject to a major taxonomic group.
    qstag = qsid.apply(tagrk)
    
    # Use majority voting to assign each query to a consensus taxon.
    mx = lambda x: max(x, key=x.tolist().count)
    QSidvote = qstag.groupby([qstag.index]).apply(mx)  # assign the Qid to the taxon that most of its Sids belong to.
    xx = QSidvote.groupby(QSidvote).count()  # summarize how many Qids of each taxon.
    xx.dropna(inplace=True)
    
    # Compile results with counts, percentages, and explode flags for visualization.
    res = pd.DataFrame(columns=['Count', 'Percentage', 'Explode'])
    res['Count'] = xx
    res.sort_values(['Count'], ascending=False, inplace=True)
    sumCount = res['Count'].sum()
    res['Percentage'] = res['Count'].div(sumCount) * 100
    res['Explode'] = [0.1 if i == target else 0 for i in res.index]
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate pie charts showing taxonomic composition of decontaminated transcripts.")
    parser.add_argument('--infile', '-i',
                        type=str,
                        help='Decontamination file ending with ".DeContX" containing DIAMOND/BLAST taxonomy data.')
    parser.add_argument('--out', '-o',
                        type=str,
                        help='Output directory for PNG and XLSX result files.')
    parser.add_argument('--db', '-d',
                        type=str,
                        help='Taxonomy database location: "remote" for online NCBI queries, or path to local taxdump directory.',
                        default="remote")
    args = parser.parse_args()
    
    # Read and preprocess the decontamination file.
    decontfile = args.infile.replace('\\', '/')
    decontheader = ['', 'sseqid', 'qstart', 'qend', 'staxids']
    Dec = pd.read_csv(decontfile, index_col=0, names=decontheader, sep="\t")
    
    # Extract genus and species from filename.
    gnm = decontfile.split('/')[-1].split('_')
    gnsp = (gnm[0].isupper() and gnm[1:4] or gnm[0:3])
    
    # Determine whether to use remote or local taxonomy database.
    if args.db.lower() == "remote":
        # Query NCBI Taxonomy online.
        ReMoTe = True
    elif all([os.path.exists(args.db + '/nodes.dmp'), os.path.exists(args.db + '/names.dmp'), os.path.exists(args.db + '/merged.dmp'), os.path.exists(args.db + '/delnodes.dmp')]):
        # Use local NCBI taxonomy dump or custom database of same format.
        ReMoTe = False
        TaxNodes = pd.read_csv(args.db + '/nodes.dmp', names=['taxid', 'parent', 'rank'], index_col=0, header=None, usecols=[0, 2, 4], sep="\t")
        TaxNames = pd.read_csv(args.db + '/names.dmp', names=['taxid', 'nametxt', 'uniquename', 'nameclass'], index_col=0, header=None, usecols=[0, 2, 4, 6], sep="\t")
        TaxNames = TaxNames[TaxNames['nameclass'] == 'scientific name']
        TaxMerge = pd.read_csv(args.db + '/merged.dmp', names=['oldid', 'validid'], index_col=0, header=None, usecols=[0, 2], sep="\t")
        TaxDel = pd.read_csv(args.db + '/delnodes.dmp', index_col=None, header=None, usecols=[0], sep="\t").iloc[:, 0]
    else:
        print('Wrong input for the argument of "--db".\n')
        exit()
    
    # Calculate taxonomic frequency distribution.
    RankRef = qclass(gnsp[0], rank='phylum', remote=ReMoTe)
    RankFreq = freq(Dec, RankRef, stepmax=50, remote=ReMoTe)
    
    # Generate and save pie chart.
    tit = RankRef + '+' + ' '.join(gnsp)
    print('Drawing ' + tit + '\n')
    taxcol = pd.Series(['darkgray', 'cyan', 'royalblue', 'brown', 'coral', 'red', 'darkgreen', 'hotpink', 'purple', 'orange', 'yellow', 'fuchsia', 'goldenrod', 'khaki', 'olive', 'orchid', 'salmon', 'wheat'], index=[RankRef, 'Arthr', 'Ecdys+', 'Spira', 'Proto+', 'Deute', 'Xenac', 'Bilat+', 'Cnida', 'Cteno', 'Placo', 'Porif', 'Virid', 'Fungi', 'Bacte', 'Archa', 'Virus', 'others'])
    col = taxcol[RankFreq.index].tolist()
    
    # Create and configure pie chart.
    plt.figure(facecolor='snow')
    plt.axes(aspect="equal")  # Ensure a perfect circle
    plt.xlim(0, 8)
    plt.ylim(0, 8)
    plt.pie(x=RankFreq['Count'], colors=col, labels=RankFreq.index, labeldistance=1.1, radius=1, explode=RankFreq['Explode'], frame=0)
    plt.title(tit.replace('+', '\n'))
    plt.xticks(())
    plt.yticks(())
    plt.savefig(args.out + '/' + decontfile.split('/')[-1].replace('.decontX', '.png'))
    
    # Generate summary statistics and save to Excel.
    RankFreq = RankFreq.append(pd.Series(RankFreq.sum(axis=0), name='Sum'))
    RankFreq['Explode'] = RankFreq['Explode'].apply(lambda x: '*' if x > 0 else None)
    RankFreq.at['Sum', 'Explode'] = None
    RankFreq.to_excel(args.out + '/' + decontfile.split('/')[-1].replace('.decontX', '.xlsx'))
