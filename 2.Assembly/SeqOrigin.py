#!/To/Your/Directory/miniconda3/bin/python3.9

"""
Generate pie charts showing taxonomic composition of assembled transcripts.

This script analyzes NCBI taxonomy for sequence hits and creates visual representations
of transcript origins by phylum and major taxonomic groups. It requires decontamination
output files with BLAST/DIAMOND taxonomy information.
"""

# Load required libraries.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import ceil
from Bio import Entrez

# Set the email address required by NCBI Entrez API.
Entrez.email = 'A.N.Other@example.com'

def qclass(GenName, rank):
    """Retrieve the taxonomic name at a specified rank for a query organism."""
    qtaxrefId = Entrez.read(Entrez.esearch(db='taxonomy', term=GenName))['IdList'][0]
    qtaxrefLin = Entrez.read(Entrez.efetch(db='taxonomy', id=qtaxrefId))[0]['LineageEx']
    qreftax = [rk['ScientificName'] for rk in qtaxrefLin if rk['Rank'].lower() == rank][0]
    return qreftax

def freq(dec, explode):
    """Analyze taxonomic composition of decontaminated sequences and generate frequency counts.
    
    Retrieves NCBI taxonomy information for all sequence hits, classifies them into
    major taxonomic groups, and performs voting-based consensus assignment.
    """
    Q = dec.index
    Qid = []
    Sid = []
    # Extract query and subject IDs from the decontamination table.
    for i in range(len(Q)):
        d = str(dec.iloc[i,]['staxids']).split(';')
        Sid.extend(d)
        Qid.extend([Q[i] for k in range(len(d))])
    # Remove null entries.
    Qid = [Qid[i] for i in range(len(Sid)) if str(Sid[i]) != 'nan']
    Sid = [Sid[i] for i in range(len(Sid)) if str(Sid[i]) != 'nan']
    QS = pd.DataFrame({'Qid': Qid, 'Sid': Sid})
    
    # Fetch taxonomy data from NCBI in batches.
    SidUni = list(set(Sid))
    stp = [x * 10000 for x in range(ceil(len(SidUni) / 10000))]
    rk = []
    for i in stp:
        handle = Entrez.efetch(db="taxonomy", retstart=i, retmax=10000, id=SidUni)
        TaxInfo = Entrez.read(handle)
        handle.close()
        rk.extend([[t['AkaTaxIds'][0] if 'AkaTaxIds' in t.keys() else t['TaxId'], t['Lineage']] for t in TaxInfo])
    
    # Remove unresolvable taxonomy IDs.
    dp = list(set(SidUni) - set([x[0] for x in rk]))
    QS = QS.loc[~np.in1d(QS['Sid'], dp)]
    QS.index = QS['Qid']

    # Classify sequences into taxonomic groups based on NCBI lineage.
    rk = np.array(rk)
    for i in range(rk.shape[0]):
        r = rk[i, 1].split('; ')
        # Assign shorthand labels for major taxonomic groups.
        if explode in r:
            rk[i, 1] = explode
        elif 'Arthropoda' in r:
            rk[i, 1] = 'Arthr'
        elif 'Ecdysozoa' in r:
            rk[i, 1] = 'Ecdys+'
        elif 'Spiralia' in r:
            rk[i, 1] = 'Spira'
        elif 'Protostomia' in r:
            rk[i, 1] = 'Proto+'
        elif 'Deuterostomia' in r:
            rk[i, 1] = 'Deute'
        elif 'Xenacoelomorpha' in r:
            rk[i, 1] = 'Xenac'
        elif 'Bilateria' in r:
            rk[i, 1] = 'Bilat+'
        elif 'Cnidaria' in r:
            rk[i, 1] = 'Cnida'
        elif 'Ctenophora' in r:
            rk[i, 1] = 'Cteno'
        elif 'Placozoa' in r:
            rk[i, 1] = 'Placo'
        elif 'Porifera' in r:
            rk[i, 1] = 'Porif'
        elif 'Viridiplantae' in r:
            rk[i, 1] = 'Virid'
        elif 'Fungi' in r:
            rk[i, 1] = 'Fungi'
        elif 'Bacteria' in r:
            rk[i, 1] = 'Bacte'
        elif 'Archaea' in r:
            rk[i, 1] = 'Archa'
        elif 'Viruses' in r:
            rk[i, 1] = 'Virus'
        else:
            rk[i, 1] = 'XXXXX'

    # Map taxonomy ranks to each sequence.
    for k in range(rk.shape[0]):
        QS.loc[np.in1d(QS['Sid'], rk[k, 0]), 'Rank'] = rk[k, 1]

    # Use majority voting to assign a consensus taxonomic group per query.
    mx = lambda x: max(x, key=x.tolist().count)
    QSidvote = QS['Rank'].groupby([QS.index]).apply(mx)
    xx = QSidvote.groupby(QSidvote).count()
    res = pd.DataFrame({'Count': xx.to_list(), 'Expl': [0.1 if i == RankRef else 0 for i in xx.index]}, index=xx.index)
    return res.sort_values(['Count'], ascending=False)

if __name__ == "__main__":
    # List of decontaminated transcriptome files to analyze.
    files = ['Armorloricus_elegans_SRS1012655_trinity.Trinity.decontX', 'Priapulus_caudatus_SRS586026_trinity.Trinity.decontX', 'Pycnophyes_sp_SRS4661771_trinity.Trinity.decontX', 'Echinoderes_dujardinii_SRS4406492_trinity.Trinity.decontX', 'Epiperipatus_trinidadensis_SRS8176931_trinity.Trinity.decontX', 'Euperipatoides_rowelli_SRS6182028_trinity.Trinity.decontX', 'Echiniscus_testudo_SRS4174635_trinity.Trinity.decontX', 'Echiniscoides_sigismundi_SRS3414361_trinity.Trinity.decontX', 'Poikilolaimus_oxycercus_SRS2523288_trinity.Trinity.decontX', 'Plectus_sambesii_SRS4078035_trinity.Trinity.decontX', 'Mesobiotus_philippinicus_SRS4118511_trinity.Trinity.decontX', 'Hypsibius_exemplaris_SRS9243152_trinity.Trinity.decontX', 'Bursaphelenchus_mucronatus_SRS3217828_trinity.Trinity.decontX', 'Cooperia_punctata_SRS4459542_trinity.Trinity.decontX', 'Trichostrongylus_colubriformis_SRS5400356_trinity.Trinity.decontX', 'Trichuris_suis_SRS2763989_trinity.Trinity.decontX', 'Poikilolaimus_oxycercus_SRS2523288_trinity.Trinity.decontX', 'Ascaridia_galli_SRS5400353_trinity.Trinity.decontX', 'Gnathostoma_spinigerum_SRS3998458_trinity.Trinity.decontX', 'Brugia_malayi_SRS5176347_trinity.Trinity.decontX', 'Aphelenchoides_fragariae_SRS3320023_trinity.Trinity.decontX', 'Halomonhystera_disjuncta_SRS748703_trinity.Trinity.decontX', 'Neocamacolaimus_parasiticus_SRS8406880_trinity.Trinity.decontX', 'Trichostrongylus_colubriformis_SRS5400356_trinity.Trinity.decontX', 'Bunonema_sp_SRS7176506_trinity.Trinity.decontX', 'Pristionchus_pacificus_SRS1608351_trinity.Trinity.decontX', 'Caenorhabditis_elegans_SRS5294341_trinity.Trinity.decontX', 'Ascaris_suum_SRS3188695_trinity.Trinity.decontX', 'Acrobeloides_nanus_SRS1818667_trinity.Trinity.decontX', 'Steinernema_feltiae_SRS7622117_trinity.Trinity.decontX', 'Meloidogyne_incognita_SRS4460030_trinity.Trinity.decontX', 'Angiostrongylus_cantonensis_SRS2803151_trinity.Trinity.decontX', 'Haemonchus_contortus_SRS2716849_trinity.Trinity.decontX', 'Pontonema_vulgare_SRS4026190_trinity.Trinity.decontX', 'Trichinella_spiralis_SRS4151078_trinity.Trinity.decontX', 'Strigamia_maritima_SRS3395075_trinity.Trinity.decontX', 'Meridionale_flava_SRS7798390_trinity.Trinity.decontX', 'Sinella_curviseta_SRS3862215_trinity.Trinity.decontX', 'Haemaphysalis_longicornis_SRS8459460_trinity.Trinity.decontX']
    
    Dir_csv = 'E:/Nema_backup/data/decontamination'
    decontheader = ['', 'sseqid', 'qstart', 'qend', 'staxids']
    
    # Define taxonomic group labels and corresponding colors for pie chart visualization.
    taxa = ['AAAAA', 'Arthr', 'Ecdys+', 'Spira', 'Proto+', 'Deute', 'Xenac', 'Bilat+', 'Cnida', 'Cteno', 'Placo', 'Porif', 'Virid', 'Fungi', 'Bacte', 'Archa', 'Virus', 'XXXXX']
    colors = ['darkgray', 'cyan', 'royalblue', 'brown', 'coral', 'red', 'darkgreen', 'hotpink', 'purple', 'orange', 'yellow', 'fuchsia', 'goldenrod', 'khaki', 'olive', 'orchid', 'salmon', 'wheat']
    
    # Process each transcriptome file.
    for fl in files:
        decontfile = Dir_csv + '/' + fl
        Dec = pd.read_csv(decontfile, index_col=0, names=decontheader, sep="\t")
        
        # Extract organism genus from filename and determine its reference phylum.
        gnm = decontfile.split('/')[-1].split('_')[0]
        RankRef = qclass(gnm, rank='phylum')
        tit = RankRef + ' ' + ' '.join(fl.split('_')[:2])
        print(tit + '\n')
        
        # Calculate taxonomic frequency distribution.
        RankFreq = freq(Dec, RankRef)
        taxa[0] = RankRef
        col = [colors[taxa.index(i)] for i in RankFreq.index]

        # Create and save pie chart.
        plt.figure(facecolor='snow')
        plt.axes(aspect="equal")  # Ensure a perfect circle
        plt.xlim(0, 8)
        plt.ylim(0, 8)
        plt.pie(x=RankFreq['Count'], colors=col, labels=RankFreq.index, labeldistance=1.1, radius=1, explode=RankFreq['Expl'], frame=0)
        plt.title(tit)
        plt.xticks(())
        plt.yticks(())
        plt.savefig('D:/Projects/BristolProjects/Nematoda/data/RNA-seq/TranscriptContents/' + tit + '.png')

