#!/To/Your/Directory/anaconda3/bin/python3.9
# PAMLphylip.py
# Helper script to convert gene alignments and partition schemes into a
# PAML-compatible relaxed phylip format for downstream MCMCTREE / PAML analyses.
#
# This script supports multiple partition-generation methods, including sequence-
# based clustering, rate-based clustering, and tree-based clustering.
# It reads an IQ-TREE-style partition file and optionally uses sequence or tree
# information to aggregate partitions, then writes the resulting alignments in
# relaxed phylip format.

import re, os
import numpy as np
from numpy.linalg import eigvals, eig, svd
import pandas as pd
# from scipy import linalg as la
from sklearn.cluster import DBSCAN, KMeans
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn_extra.cluster import KMedoids
from kneed import KneeLocator
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

# Parse an IQ-TREE/NEXUS charset block and return partition start/end indices.
# The input lines should look like: "charset gene1 = 1-100".
# The returned DataFrame uses zero-based coordinates for Python slicing.
def GetPartition(old_partitions):
    CharSet = pd.DataFrame([line.strip('charset ').replace('  ', ' ').split(' = ') for line in old_partitions if line.startswith('charset')])
    CharSet.columns = ['charset', 'partition']
    CharSet.set_index('charset', inplace=True)
    CharSet[['start','end']] = CharSet['partition'].apply(lambda x: pd.Series([int(i) for i in x.split('-')]))
    CharSet['start'] = CharSet['start']-1
    return CharSet[['start','end']]

# Merge partition definitions by substitution model, grouping charsets that share the same model.
# This is useful when the input partition file contains both charset definitions and model assignments.
def GetPartitionModel(old_partitions):
    CharSet = pd.DataFrame([line.strip('charset ').replace('  ', ' ').split(' = ') for line in old_partitions if line.startswith('charset')])
    CharSet.columns = ['charset', 'partition']
    CharSet.set_index('charset', inplace=True)
    Model = pd.DataFrame([line.split(': ') for line in old_partitions if ": " in line])
    Model.columns = ['model','charset']
    Model.set_index('charset', inplace=True)
    X = pd.concat([CharSet, Model], axis=1)
    MergeX = X.groupby('model').apply(lambda x: list(map(int, re.split('[ -]', ' '.join(x['partition'].tolist())))))
    DataSet = pd.DataFrame(MergeX.apply(lambda x: [x[0::2], x[1::2]]).to_dict()).transpose()
    DataSet.columns = ['start', 'end']
    DataSet['start'] = DataSet['start'].apply(lambda x: [i-1 for i in x])
    return DataSet

# Compute the von Neumann entropy for a sequence of discrete states.
# This is used to measure sequence variability over a set of positions.
def von_neumann_entropy(x, pos):
    x = np.array(x).reshape(-1)
    x = x[np.where(np.isin(x,pos))]
    totlen = len(x)
    freq = np.zeros(len(pos))
    for i in range(len(pos)):
        freq[i] = np.sum(x == pos[i])/totlen
    rho = np.matrix(np.diag(freq)) * np.diag([1]*len(pos))
    mu = 1 / np.trace(rho)
    rho_norm = mu * rho
    diag = eigvals(rho_norm) * np.log(eigvals(rho_norm)) / np.log(20) 
    return -np.sum(diag[~np.isnan(diag)])

# Calculate entropy of an amino acid column using a substitution matrix.
# Returns weighted entropy and the proportion of nongap residues.
def von_neumann_entropy_prot(char, sub_mat='BLOSUM62'):
    substmatx = np.matrix(substitution_matrices.load(sub_mat))[:20,:20]
    pos = np.array(list(substitution_matrices.load(sub_mat).alphabet[:20]))
    char = np.array(list(char)).reshape(-1)
    p_nongap = np.sum(np.isin(char,pos))/len(char)
    char = char[np.where(np.isin(char,pos))]
    totlen = len(char)
    freq = np.zeros(len(pos))
    for i in range(len(pos)):
        freq[i] = np.sum(char == pos[i])/totlen
    rho = np.matrix(np.diag(freq)).dot(substmatx)
    # rho = np.matrix(np.diag(freq) * substmatx) # this makes no difference to the result than the above line
    mu = 1 / np.trace(rho)
    rho_norm = mu * rho
    # diag = np.diag(rho_norm.dot(la.logm(rho_norm) / np.log(20))) # this is the same as the below line
    diag = eigvals(rho_norm) * np.log(eigvals(rho_norm)) / np.log(20) 
    return -np.sum(diag[~np.isnan(diag)])*p_nongap, p_nongap

# Apply a sliding window across protein positions and compute mean entropy.
def von_neumann_entropy_prot_win(Seq, w, sub_mat='BLOSUM62'):
    w = int(w)
    seq_len = len(Seq[0])
    slice = zip([0] * w + list(range(seq_len - w)), list(range(seq_len)[w+1:]) + [seq_len] * (w + 1))
    vne = []
    for start, end in slice:
        hc = np.array([von_neumann_entropy_prot(Seq[:, i], sub_mat) for i in range(start, end)])
        x = np.sum(hc, axis=0)
        vne.append(x[0]/x[1])
    return np.array(vne)

# Compute average sequence and site entropy across a gene alignment.
def hv_vne_mean(Seq, sub_mat='BLOSUM62'):
    Seq = MultipleSeqAlignment([i.upper() for i in Seq if not set(i.seq).issubset(set(['X','Z','B','-','?','*']))])
    h_vne = pd.Series([Seq[i].seq for i in range(len(Seq))]).apply(lambda x: von_neumann_entropy_prot(x, sub_mat))
    h_vne_mean = h_vne.apply(lambda x: x[0]/x[1]).mean()
    v_vne = pd.Series([Seq[:, i] for i in range(len(Seq[0]))]).apply(lambda x: von_neumann_entropy_prot(x, sub_mat))
    v_vne_mean = v_vne.apply(lambda x: x[0]/x[1]).mean()
    return pd.Series([h_vne_mean, v_vne_mean], index=['h_vne_mean', 'v_vne_mean'])

# Convert a rate table into a per-partition summary for clustering.
def gene_rate_2x(rate_table):
    rate_median = rate_table['Rate'].groupby(rate_table.index).median()
    all_cat = rate_table['Cat'].drop_duplicates()
    rate_cat_vne = rate_table['Cat'].groupby(rate_table.index).apply(lambda x: von_neumann_entropy(x, all_cat))
    median_vne = pd.concat([rate_median, rate_cat_vne], axis=1)
    median_vne.columns = ['rate_median', 'rate_cat_vne']    
    return median_vne

# Calculate the log-median branch length for each gene tree.
def tree_BL_median(trees):
    tree_BL_median = list(map(lambda x: np.log(np.median([i.branch_length for i in x.get_terminals()])), trees))
    return np.array(tree_BL_median)

# Compute a distance matrix for a tree based on taxon patristic distances.
def tree_dist(tree, taxa):
    tx = np.array([i.name for i in tree.get_terminals()])
    ntaxa = len(taxa)
    tree_dist = list(map(lambda x: np.array([tree.distance(x, i) if i in tx else 0 for i in taxa] if x in tx else [0] * ntaxa), taxa))
    return np.array(tree_dist)

# Convert tree distance information into low-dimensional SVD features for clustering.
def tree_dist_svg(trees, taxa, genes, nsv=2):
    tr_dist_eig = np.ones((len(trees), len(taxa)))
    for i in range(len(trees)):
        tr_dist = np.sqrt(tree_dist(trees[i], taxa))
        val, vec = eig(tr_dist)
        tr_dist_eig[i,:] = vec[:, 0]
    U, D, VT = svd(tr_dist_eig)
    svv = pd.DataFrame(U[:, :nsv], index=genes)
    return svv

# Automatic cluster selection using KMeans and silhouette score.
def kmeans(data, max_clusters = 20):
    k_series = pd.Series(range(2, max_clusters + 1), index=range(2, max_clusters + 1))
    kmeans_fit_series = k_series.apply(lambda x: KMeans(n_clusters=x, random_state=50, init="k-means++").fit(data))
    silhouette_series = kmeans_fit_series.apply(lambda x: silhouette_score(data, x.labels_, metric="euclidean",sample_size=1000,random_state=200))
    # print(silhouette_series)
    return kmeans_fit_series[silhouette_series.idxmax()]

# Automatic cluster selection using KMedoids and silhouette score.
def kmedoids(data, max_clusters = 20):
    k_series = pd.Series(range(2, max_clusters + 1), index=range(2, max_clusters + 1))
    kmedoids_fit_series = k_series.apply(lambda x: KMedoids(n_clusters=x, random_state=50, init="heuristic", method='pam').fit(data))
    silhouette_series = kmedoids_fit_series.apply(lambda x: silhouette_score(data, x.labels_, metric="euclidean",sample_size=1000,random_state=200))
    # print(silhouette_series)
    return kmedoids_fit_series[silhouette_series.idxmax()]

# DBSCAN clustering with optional automatic parameter tuning.
def dbscan(data, eps=None, min_samples=None):
    # compute minimum number of samples in a cluster when not specified by user
    if min_samples == None:
        min_samples = 2 * np.ndim(data)
    # compute optimized epsilon when not specified by user   
    if eps == None: 
        neighbors = NearestNeighbors(n_neighbors=20, algorithm='ball_tree').fit(data)
        distances, indices = neighbors.kneighbors(data)
        distances = np.sort(distances[:,10], axis=0)
        knee = KneeLocator(np.arange(len(distances)), distances, S=1, curve='convex', direction='increasing', interp_method='polynomial')
        eps = distances[knee.knee]
    # compute DBSCAN with optimized epsilon and minimum number of samples
    db_cluster = DBSCAN(eps=eps, min_samples=min_samples).fit(data)
    # len(set(db_cluster.labels_))-(1 if -1 in db_cluster.labels_ else 0)
    # v_measure_score(y, db_cluster.labels_)
    return db_cluster

# Rebuild partition definitions using cluster labels and original loci coordinates.
def PartitionTransformer(loci, loci_loc, clustering_model):
    grp = pd.DataFrame({'loci':loci, 'start':loci_loc['start'][loci].values, 'end':loci_loc['end'][loci].values}, index=clustering_model.labels_).sort_index()
    NewPart = grp.groupby(grp.index).apply(lambda x: pd.Series({'start':list(x['start']),'end':list(x['end'])}))
    return NewPart

# Concatenate a list of alignment slices into a single combined alignment.
def AlignAddition(seq):
    CombAlign = MultipleSeqAlignment([SeqRecord(Seq(''), id=i.id) for i in seq[0]])
    for i in seq:
        CombAlign += i
    return CombAlign

# Create new alignment partitions from original alignments and partition coordinates.
def SepAlign(old_align, new_partitions):
    new_align = new_partitions.apply(lambda x: AlignAddition([old_align[:, st:ed] for st, ed in zip(x['start'],x['end'])]), axis=1)
    return new_align

if __name__ == "__main__":
    # Parse command-line arguments for input alignments, partition definitions,
    # auxiliary tree/rate data, clustering method, and output destination.
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq', '-s',
                        type=str,
                        help='Input sequence file (multiple alignments) in phylip format.')
    parser.add_argument('--partition', '-p',
                        type=str,
                        help='Partition file in IQTREE nexus format.',
                        default='')
    parser.add_argument('--infd', '-i',
                        type=str,
                        help='Input multi-tree file in newick format or a directory containing all gene trees when using tree-based clustering. It can also be a tab-delimited IQTREE site-rate output for rate-based clustering. Leave empty for seq-based methods.',
                        default='')
    parser.add_argument('--method', '-m',
                        type=str,
                        help='Partitioning method. If omitted, the script uses the provided IQ-TREE partition file and merges partitions by shared models. Supported methods: seq-based-kmeans, seq-based-kmedoids, seq-based-dbscan, rate-based-kmeans, rate-based-kmedoids, rate-based-dbscan, tree-based-kmeans, tree-based-kmedoids, tree-based-dbscan, tree-based-svd-kmeans, tree-based-svd-kmedoids, tree-based-svd-dbscan.',
                        default='')
    parser.add_argument('--out', '-o',
                        type=str,
                        help='Output file in PAML phylip-relaxed format.')
    args = parser.parse_args()

    # Normalize file paths for cross-platform compatibility.
    InputAlignFile = args.seq.replace('\\','/')
    PartitionFile = args.partition.replace('\\','/')
    InputExtraFile = args.infd.replace('\\','/')
    OutputFile = args.out.replace('\\','/')

    # Read the relaxed phylip alignment and extract taxon names.
    PhylipAlignments = AlignIO.read(InputAlignFile, "phylip-relaxed")
    taxa = np.array([i.id for i in PhylipAlignments])

    if os.path.isfile(PartitionFile):
        with open(PartitionFile, 'r') as PartSch:
            PartitionScheme = [line.strip('\n').strip(',').strip(';').strip(' ') for line in PartSch.readlines()]
        PartSch.close()
        del PartSch

        # Each method constructs a new partition scheme using different criteria.
        if args.method.lower() == 'seq-based-kmeans':
            GeneLoc = GetPartition(PartitionScheme)
            gene_vne = GeneLoc.apply(lambda x: hv_vne_mean(PhylipAlignments[:, x[0]:x[1]], sub_mat='BLOSUM90'), axis=1)
            seq_kmeans_model = kmeans(gene_vne, max_clusters=10)
            NewPart = PartitionTransformer(gene_vne.index, GeneLoc, seq_kmeans_model)
        if args.method.lower() == 'seq-based-kmedoids':
            GeneLoc = GetPartition(PartitionScheme)
            gene_vne = GeneLoc.apply(lambda x: hv_vne_mean(PhylipAlignments[:, x[0]:x[1]], sub_mat='BLOSUM90'), axis=1)
            seq_kmedoids_model = kmedoids(gene_vne, max_clusters=10)
            NewPart = PartitionTransformer(gene_vne.index, GeneLoc, seq_kmedoids_model)        
        elif args.method.lower() == 'seq-based-dbscan':
            GeneLoc = GetPartition(PartitionScheme)
            gene_vne = GeneLoc.apply(lambda x: hv_vne_mean(PhylipAlignments[:, x[0]:x[1]], sub_mat='BLOSUM90'), axis=1)
            seq_dbscan_model = dbscan(gene_vne)
            NewPart = PartitionTransformer(gene_vne.index, GeneLoc, seq_dbscan_model)
        elif args.method.lower() == 'rate-based-kmeans':
            GeneLoc = GetPartition(PartitionScheme)
            with open(InputExtraFile, 'r') as RateFile:
                RateTable = [line.strip('\n').strip(',').strip(';').strip(' ').split('\t') for line in RateFile.readlines() if not line.startswith('#')]
            RateFile.close()
            del RateFile
            site_rate = list(map(lambda x: [GeneLoc.index[int(x[0])-1], int(x[1]), float(x[2]), int(x[3]), float(x[4])], RateTable[1:]))
            site_rate = pd.DataFrame(site_rate, columns=RateTable[0])
            site_rate.set_index(site_rate['Part'], inplace=True)
            site_rate = site_rate.iloc[:, 1:]
            rate_med_catvne = gene_rate_2x(site_rate)
            rate_kmeans_model = kmeans(rate_med_catvne, max_clusters=10)
            NewPart = PartitionTransformer(rate_med_catvne.index, GeneLoc, rate_kmeans_model)        
        elif args.method.lower() == 'rate-based-kmedoids':
            GeneLoc = GetPartition(PartitionScheme)
            with open(InputExtraFile, 'r') as RateFile:
                RateTable = [line.strip('\n').strip(',').strip(';').strip(' ').split('\t') for line in RateFile.readlines() if not line.startswith('#')]
            RateFile.close()
            del RateFile
            site_rate = list(map(lambda x: [GeneLoc.index[int(x[0])-1], int(x[1]), float(x[2]), int(x[3]), float(x[4])], RateTable[1:]))
            site_rate = pd.DataFrame(site_rate, columns=RateTable[0])
            site_rate.set_index(site_rate['Part'], inplace=True)
            site_rate = site_rate.iloc[:, 1:]
            rate_med_catvne = gene_rate_2x(site_rate)
            rate_kmedoids_model = kmedoids(rate_med_catvne, max_clusters=10)
            NewPart = PartitionTransformer(rate_med_catvne.index, GeneLoc, rate_kmedoids_model)        
        elif args.method.lower() == 'rate-based-dbscan':
            GeneLoc = GetPartition(PartitionScheme)
            with open(InputExtraFile, 'r') as RateFile:
                RateTable = [line.strip('\n').strip(',').strip(';').strip(' ').split('\t') for line in RateFile.readlines() if not line.startswith('#')]
            RateFile.close()
            del RateFile
            site_rate = list(map(lambda x: [GeneLoc.index[int(x[0])-1], int(x[1]), float(x[2]), int(x[3]), float(x[4])], RateTable[1:]))
            site_rate = pd.DataFrame(site_rate, columns=RateTable[0])
            site_rate.set_index(site_rate['Part'], inplace=True)
            site_rate = site_rate.iloc[:, 1:]
            rate_med_catvne = gene_rate_2x(site_rate)
            rate_med_catvne_dbscan_model = dbscan(rate_med_catvne)
            NewPart = PartitionTransformer(rate_med_catvne.index, GeneLoc, rate_med_catvne_dbscan_model)
        elif args.method.lower() == 'tree-based-kmeans':
            GeneLoc = GetPartition(PartitionScheme)
            GeneTrees=list(Phylo.parse(InputExtraFile, "newick"))
            if all([sum([1 for i in PhylipAlignments[:,GeneLoc.iloc[t,0]:GeneLoc.iloc[t,1]] if not set(i.seq)==set(['?'])]) == len(GeneTrees[t].get_terminals()) for t in range(len(GeneTrees))]):
                pass
            else:
                print('The gene trees should be in same order as genes appear in the partition file. Please check the input tree file.')
                exit()
            gene_tree_median = pd.DataFrame({'x':np.ones(len(GeneTrees)), 'BL_median':tree_BL_median(GeneTrees)}, index=GeneLoc.index)
            gene_tree_median_kmeans_model = kmeans(gene_tree_median, max_clusters=10)
            NewPart = PartitionTransformer(gene_tree_median.index, GeneLoc, gene_tree_median_kmeans_model)        
        elif args.method.lower() == 'tree-based-kmedoids':
            GeneLoc = GetPartition(PartitionScheme)
            GeneTrees=list(Phylo.parse(InputExtraFile, "newick"))
            if all([sum([1 for i in PhylipAlignments[:,GeneLoc.iloc[t,0]:GeneLoc.iloc[t,1]] if not set(i.seq)==set(['?'])]) == len(GeneTrees[t].get_terminals()) for t in range(len(GeneTrees))]):
                pass
            else:
                print('The gene trees should be in same order as genes appear in the partition file. Please check the input tree file.')
                exit()
            gene_tree_median = pd.DataFrame({'x':np.ones(len(GeneTrees)), 'BL_median':tree_BL_median(GeneTrees)}, index=GeneLoc.index)
            gene_tree_median_kmedoids_model = kmedoids(gene_tree_median, max_clusters=10)
            NewPart = PartitionTransformer(gene_tree_median.index, GeneLoc, gene_tree_median_kmedoids_model)        
        elif args.method.lower() == 'tree-based-dbscan':
            GeneLoc = GetPartition(PartitionScheme)
            GeneTrees=list(Phylo.parse(InputExtraFile, "newick"))
            if all([sum([1 for i in PhylipAlignments[:,GeneLoc.iloc[t,0]:GeneLoc.iloc[t,1]] if not set(i.seq)==set(['?'])]) == len(GeneTrees[t].get_terminals()) for t in range(len(GeneTrees))]):
                pass
            else:
                print('The gene trees should be in same order as genes appear in the partition file. Please check the input tree file.')
                exit()
            gene_tree_median = pd.DataFrame({'x':np.ones(len(GeneTrees)), 'BL_median':tree_BL_median(GeneTrees)}, index=GeneLoc.index)
            gene_tree_median_dbscan_model = dbscan(gene_tree_median)
            NewPart = PartitionTransformer(gene_tree_median.index, GeneLoc, gene_tree_median_dbscan_model)        
        elif args.method.lower() == 'tree-based-svd-kmeans':
            GeneLoc = GetPartition(PartitionScheme)
            GeneTrees=list(Phylo.parse(InputExtraFile, "newick"))
            if all([sum([1 for i in PhylipAlignments[:,GeneLoc.iloc[t,0]:GeneLoc.iloc[t,1]] if not set(i.seq)==set(['?'])]) == len(GeneTrees[t].get_terminals()) for t in range(len(GeneTrees))]):
                pass
            else:
                print('The gene trees should be in same order as genes appear in the partition file. Please check the input tree file.')
                exit()
            svv = tree_dist_svg(GeneTrees, taxa, GeneLoc.index, nsv=2)
            svv_kmeans_model = kmeans(svv, max_clusters=10)
            NewPart = PartitionTransformer(svv.index, GeneLoc, svv_kmeans_model)
        elif args.method.lower() == 'tree-based-svd-kmedoids':
            GeneLoc = GetPartition(PartitionScheme)
            GeneTrees=list(Phylo.parse(InputExtraFile, "newick"))
            if all([sum([1 for i in PhylipAlignments[:,GeneLoc.iloc[t,0]:GeneLoc.iloc[t,1]] if not set(i.seq)==set(['?'])]) == len(GeneTrees[t].get_terminals()) for t in range(len(GeneTrees))]):
                pass
            else:
                print('The gene trees should be in same order as genes appear in the partition file. Please check the input tree file.')
                exit()
            svv = tree_dist_svg(GeneTrees, taxa, GeneLoc.index, nsv=2)
            svv_kmedoids_model = kmedoids(svv, max_clusters=10)
            NewPart = PartitionTransformer(svv.index, GeneLoc, svv_kmedoids_model)            
        elif args.method.lower() == 'tree-based-svd-dbscan':
            GeneLoc = GetPartition(PartitionScheme)
            GeneTrees=list(Phylo.parse(InputExtraFile, "newick"))
            if all([sum([1 for i in PhylipAlignments[:,GeneLoc.iloc[t,0]:GeneLoc.iloc[t,1]] if not set(i.seq)==set(['?'])]) == len(GeneTrees[t].get_terminals()) for t in range(len(GeneTrees))]):
                pass
            else:
                print('The gene trees should be in same order as genes appear in the partition file. Please check the input tree file.')
                exit()
            svv = tree_dist_svg(GeneTrees, taxa, GeneLoc.index, nsv=3)
            svv_dbscan_model = dbscan(svv)
            NewPart = PartitionTransformer(svv.index, GeneLoc, svv_dbscan_model)            
        else:
            NewPart = GetPartitionModel(PartitionScheme)
        NewAlignments = SepAlign(PhylipAlignments, NewPart)
    else:
        # If no partition file is provided, preserve the original alignment as a single partition.
        NewAlignments = [MultipleSeqAlignment([SeqRecord(i.seq, id=i.id) for i in PhylipAlignments])]

    # Interleaved formats (not recommended, as it may not be compatible with PAML)
    # AlignIO.write(NewAlignments.to_list(), "OutputFile", "phylip-relaxed")
    
    # Write the resulting partitions in relaxed phylip sequential format.
    with open(OutputFile, 'w') as OutFile:
        for part in NewAlignments:
            OutFile.write(' ' + str(len(part)) + ' ' + str(len(part[0])) + '\n')
            OutFile.write('\n'.join([str(i.id)+'\n'+str(i.seq) for i in part]) + '\n\n')
    OutFile.close()
    print('Converting to Phylip-format for MCMCTREE use is DONE!')
