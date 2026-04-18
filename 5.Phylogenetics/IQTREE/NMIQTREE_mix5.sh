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
IQTREE=/sw/apps/iqtree-2.1.3-Linux/bin/iqtree2
# IQTREE=/public4/home/sc56777/apps/iqtree2/bin/iqtree2

DIR_phylo=/public4/home/sc56777/Data/nmtd/TBDL.nom.rmsp
# OutGroup=PF_Amphimedon_queenslandica

cd $DIR_phylo
DS=TBDL.nom.rmsp

# homogeneous model: LG+F+R guide tree
if [[ ! -f $DIR_phylo/Nema.$DS.LG.iqtree ]]; then 
    (srun --job-name=NMIQTmix --exclusive --partition=amd_256 --nodes=1 --ntasks-per-node=1 --cpus-per-task=$NP $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+F+R --rate --mlrate -B 1000 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --runs $nruns --prefix $DIR_phylo/Nema.$DS.LG > $DIR_phylo/Nema.$DS.LG.full.log 2>&1 &) 
fi

# mixture model: LG4X guide tree
if [[ ! -f $DIR_phylo/Nema.$DS.LG4X.iqtree ]]; then 
    (srun --job-name=NMIQTmix --exclusive --partition=amd_256 --nodes=1 --ntasks-per-node=1 --cpus-per-task=$NP $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG4X+FO --rate --mlrate -B 1000 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --runs $nruns --prefix $DIR_phylo/Nema.$DS.LG4X > $DIR_phylo/Nema.$DS.LG4X.full.log 2>&1 &) 
fi

# mixture model: C60 mixture site distribution profile
if [[ ! -f $DIR_phylo/Nema.$DS.C60.iqtree ]]; then 
    (srun --job-name=NMIQTmix --exclusive --partition=amd_2T --nodes=1 --ntasks-per-node=1 --cpus-per-task=$NP $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m C60+FO+R --rate --mlrate -B 1000 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --runs 1 --prefix $DIR_phylo/Nema.$DS.C60 >> $DIR_phylo/Nema.$DS.C60.full.log 2>&1 &)
fi 
    
# mixture model: LGC60 PMSF mixture site distribution profile
[[ -f $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof.runs ]] && r0=$[$(awk 'END {print $1}' $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof.runs)+1] || r0=1 
        
if [[ $r0 -le $nruns ]]; then 
    for run in `seq $r0 $nruns` 
    do 
        $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+C60+FO+R -ft $DIR_phylo/Nema.$DS.LG4X.treefile -n 0 -T $NP --prefix $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof$run >> $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof.full.log 2>&1 && 
        curlik=$(grep 'Log-likelihood of the tree: ' Nema.$DS.LGC60.PMSF.siteprof$run.iqtree | awk '{print $5}') && 
        echo -e "${run} ${curlik}" >> $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof.runs 
    done
fi 
    
# mixture model: LGC60 PMSF tree and bootstrap
maxrow=$(awk 'BEGIN{max = ""}{if($2 > max) {max = $2; id = $1}}END{print id}' $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof.runs) 

if [[ ! -f $DIR_phylo/Nema.$DS.LGC60.PMSF.B1000.iqtree ]]; then 
    $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+C60+FO+R -fs $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof${maxrow}.sitefreq -B 1000 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --prefix $DIR_phylo/Nema.$DS.LGC60.PMSF.B1000 > $DIR_phylo/Nema.$DS.LGC60.PMSF.B1000.full.log 2>&1 
fi 

if [[ ! -f $DIR_phylo/Nema.$DS.LGC60.PMSF.b100.iqtree ]]; then 
    $IQTREE -s $DIR_phylo/FullData184taxa171097AA.phy --seqtype AA -m LG+C60+FO+R -fs $DIR_phylo/Nema.$DS.LGC60.PMSF.siteprof${maxrow}.sitefreq -b 100 --wbtl --alrt 1000 --abayes --lbp 1000 -T $NP --prefix $DIR_phylo/Nema.$DS.LGC60.PMSF.b100 > $DIR_phylo/Nema.$DS.LGC60.PMSF.b100.full.log 2>&1 
fi 	
