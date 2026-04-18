#!/bin/bash -l
#SBATCH --job-name=NMMetaEuk
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128GB
#SBATCH --exclusive

NP=64

DIR_cnd=/public1/home/scb4616/apps/anaconda3

source $DIR_cnd/bin/activate $DIR_cnd/envs/makerenv

if [ ! -d /public1/home/scb4616/Data/nmtd/preds/metaeuk ]; then
    mkdir /public1/home/scb4616/Data/nmtd/preds/metaeuk
fi

cd /public1/home/scb4616/Data/nmtd/preds/metaeuk

metaeuk easy-predict /public1/home/scb4616/Data/nmtd/assembly/Bunonema_sp_GCA_014805405.1.fna /public1/home/scb4616/db/OrthoDB/refseq_db.nematoda.metaeukdb  Bunonema_sp_GCA_014805405.1.metaeuk /public1/home/scb4616/tmp --threads $NP --remove-tmp-files 1 --num-iterations 5  

# metaeuk easy-predict /To/Your/Directory/nematoda/data/assembly/Bunonema_sp_GCA_014805405.1.fna /To/Your/Directory/db/UniRef50/uniref50.metaeukdb Bunonema_sp_GCA_014805405.1.metaeuk /To/Your/Directory/tmp --threads $NP --remove-tmp-files 1 --protein 1 --num-iterations 5 --exhaustive-search 1 

