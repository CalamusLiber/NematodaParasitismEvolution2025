#!/bin/bash
#SBATCH --job-name=NMIQTREE
#SBATCH --exclusive
#SBATCH --partition=mwvdk
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=350GB

NP=$SLURM_CPUS_PER_TASK # 32
NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE] # 4

module load apps/iqtree/2.1.3
IQTREE=/sw/apps/iqtree-2.1.3-Linux/bin/iqtree2

DIR_matrix=$DIR_msa/matrix
DIR_phylo=/To/Your/Directory/nematoda/phylo
OutGroup=PF_Amphimedon_queenslandica

cd $DIR_phylo

echo -e "#######################\nBuilding Maximum Likelihood Phylogenetic Trees Using IQTREE2.1.3...\nDataset: \n$(date)\n#######################\n\n" >> $DIR_phylo/TreeBuilding.log && 
# no partition
if [[ ! -f $DIR_phylo/Nema.$DS.np.contree ]]; then 
    echo -e "=======================\nNo Partition start...\n$(date)\n=======================\n\n" >> $DIR_phylo/TreeBuilding.log && 

    $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m MFP --msub nuclear --mset LG,WAG --rate --mlrate -B 1000 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --runs 4 -o $OutGroup --prefix $DIR_phylo/Nema.np && 

    echo -e "=======================\nNo Partition finished.\n$(date)\n=======================\n\n" >> $DIR_phylo/TreeBuilding.log 
fi



