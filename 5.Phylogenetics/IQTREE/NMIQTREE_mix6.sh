#!/bin/bash
#SBATCH --job-name=NMIQTmix
#SBATCH --exclusive
#SBATCH --partition=amd_2T
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=64
#SBATCH --mem=360GB

NP=64 # 16
nruns=4
# NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE] # 4

# module load apps/iqtree/2.1.3
# IQTREE=/sw/apps/iqtree-2.1.3-Linux/bin/iqtree2
# IQTREE=/public4/home/sc56777/apps/iqtree2/bin/iqtree2
# IQTREE=/apps/iqtree2/bin/iqtree2
IQTREE=/mnt/data0/lvliang/iqtree2/bin/iqtree2

# DIR_phylo=/public4/home/sc56777/Data/nmtd/TBDL.nom.rmsp.Mixed.Tree
# DIR_phylo=/mnt/data/TBDL.nom.rmsp.Mixed.Tree
DIR_phylo=/To/Your/Directory/nematoda/phylo/TBDL.nom.rmsp.Mixed.Tree
# OutGroup=PF_Amphimedon_queenslandica

cd $DIR_phylo
DS=Mixed.Tree

# homogeneous model: LG+F+R partially constrained tree (as PMSF guide tree)
if [[ ! -f $DIR_phylo/Nema.$DS.LG.constr.iqtree ]]; then 
    $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+F+R -g $DIR_phylo/AllBackbone.bp.merged.constraint.tre --rate --mlrate -B 1000 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --runs $nruns --prefix $DIR_phylo/Nema.$DS.LG.constr > $DIR_phylo/Nema.$DS.LG.constr.full.log 2>&1 && 
    cp $DIR_phylo/Nema.$DS.LG.constr.treefile $DIR_phylo/topology.test 
fi

# mixture model: cLG PMSF mixture site distribution profile (run on BSCC)
[[ -f $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof.runs ]] && r0=$[$(awk 'END {print $1}' $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof.runs)+1] || r0=1 
        
if [[ $r0 -le $nruns ]]; then 
    for run in `seq $r0 $nruns` 
    do 
        $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+C60+FO+R -ft $DIR_phylo/Nema.$DS.LG.constr.treefile -n 0 -T $NP --prefix $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof$run >> $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof.full.log 2>&1 && 
        curlik=$(grep 'Log-likelihood of the tree: ' Nema.$DS.cLG.PMSF.siteprof$run.iqtree | awk '{print $5}') && 
        echo -e "${run} ${curlik}" >> $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof.runs 
    done
fi 
    
# mixture model: cLG PMSF tree and bootstrap
maxrow=$(awk 'BEGIN{max = ""}{if($2 > max) {max = $2; id = $1}}END{print id}' $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof.runs) 

if [[ ! -f $DIR_phylo/Nema.$DS.cLG.PMSF.B1000.iqtree ]]; then 
    $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+C60+FO+R -fs $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof${maxrow}.sitefreq -B 1000 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --prefix $DIR_phylo/Nema.$DS.cLG.PMSF.B1000 > $DIR_phylo/Nema.$DS.cLG.PMSF.B1000.full.log 2>&1 && 
    cp $DIR_phylo/Nema.$DS.cLG.PMSF.B1000.treefile $DIR_phylo/topology.test 
fi 

if [[ ! -f $DIR_phylo/Nema.$DS.cLG.PMSF.b100.iqtree ]]; then 
    $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+C60+FO+R -fs $DIR_phylo/Nema.$DS.cLG.PMSF.siteprof${maxrow}.sitefreq -b 100 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --prefix $DIR_phylo/Nema.$DS.cLG.PMSF.b100 > $DIR_phylo/Nema.$DS.cLG.PMSF.b100.full.log 2>&1 && 
    cp $DIR_phylo/Nema.$DS.cLG.PMSF.b100.treefile $DIR_phylo/topology.test 
fi 	

# tree/topology test
cd $DIR_phylo/topology.test 

# multiple test
# Nema.TBDL.nom.rmsp.C60.PMSF.b100.treefile Nema.Mixed.Tree.cLG.PMSF.b100.treefile 
cat Nema.TBDL.nom.rmsp.LG.treefile Nema.TBDL.nom.rmsp.LG4X.treefile Nema.TBDL.nom.rmsp.C60.treefile Nema.TBDL.nom.rmsp.C60.PMSF.B1000.treefile Nema.Mixed.Tree.cLG.treefile Nema.Mixed.Tree.cLG.PMSF.B1000.treefile > Nema.All.Trees.treelist && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+F+R -z $DIR_phylo/topology.test/Nema.All.Trees.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.All.Trees.treetest1 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+C60+FO+R -z $DIR_phylo/topology.test/Nema.All.Trees.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.All.Trees.treetest2 

# pairwise test
# LG+F+R
cat Nema.TBDL.nom.rmsp.LG.treefile Nema.Mixed.Tree.cLG.treefile > Nema.LG-cLG.treelist && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+F+R -z $DIR_phylo/topology.test/Nema.LG-cLG.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.LG-cLG.treetest1 && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+C60+FO+R -z $DIR_phylo/topology.test/Nema.LG-cLG.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.LG-cLG.treetest2 

cat Nema.TBDL.nom.rmsp.C60.PMSF.B1000.treefile Nema.Mixed.Tree.cLG.PMSF.B1000.treefile > Nema.C60-cLG.PMSF.B1000.treelist && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+F+R -z $DIR_phylo/topology.test/Nema.C60-cLG.PMSF.B1000.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.C60-cLG.PMSF.B1000.treetest1 && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+C60+FO+R -z $DIR_phylo/topology.test/Nema.C60-cLG.PMSF.B1000.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.C60-cLG.PMSF.B1000.treetest2 

cat Nema.TBDL.nom.rmsp.C60.PMSF.B1000.treefile Nema.Mixed.Tree.cLG.treefile > Nema.C60.PMSF-cLG.treelist && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+F+R -z $DIR_phylo/topology.test/Nema.C60.PMSF-cLG.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.C60.PMSF-cLG.treetest1 && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+C60+FO+R -z $DIR_phylo/topology.test/Nema.C60.PMSF-cLG.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.C60.PMSF-cLG.treetest2 

cat Nema.Mixed.Tree.cLG.treefile Nema.TBDL.nom.rmsp.cLG.PMSF.B1000.treefile > Nema.cLG-cLG.PMSF.treelist && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+F+R -z $DIR_phylo/topology.test/Nema.cLG-cLG.PMSF.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.cLG-cLG.PMSF.treetest1 && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+C60+FO+R -z $DIR_phylo/topology.test/Nema.cLG-cLG.PMSF.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.cLG-cLG.PMSF.treetest2 

cat Nema.TBDL.nom.rmsp.LG.treefile Nema.TBDL.nom.rmsp.C60.PMSF.B1000.treefile > Nema.LG-C60.PMSF.treelist && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+F+R -z $DIR_phylo/topology.test/Nema.LG-C60.PMSF.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.LG.C60.PMSF.treetest1 && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+C60+FO+R -z $DIR_phylo/topology.test/Nema.LG-C60.PMSF.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.LG.C60.PMSF.treetest2 

cat Nema.TBDL.nom.rmsp.LG.treefile Nema.TBDL.nom.rmsp.cLG.PMSF.B1000.treefile > Nema.LG-cLG.PMSF.treelist && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+F+R -z $DIR_phylo/topology.test/Nema.LG-cLG.PMSF.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.LG-cLG.PMSF.treetest1 && 
$IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy -m LG+C60+FO+R -z $DIR_phylo/topology.test/Nema.LG-cLG.PMSF.treelist -n 0 -zb 100000 -zw -au -T $NP --prefix $DIR_phylo/topology.test/Nema.LG-cLG.PMSF.treetest2 



