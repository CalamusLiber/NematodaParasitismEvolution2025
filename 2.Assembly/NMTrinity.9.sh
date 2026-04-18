#!/bin/sh

cd /To/Your/Directory/nematoda/data
DIR_cur=$(pwd)
chmod -R 777 $DIR_cur/NMTrinity.9.x.sh

rnaseq_list=$(ls ./rnaseq/*.R1.fastq | sed -e "s#.\/rnaseq\/##g;s#.R1.fastq##g")

for nrna in $rnaseq_list
do
    
    sbatch --partition=mwvdk --nodes=1 --ntasks=32 --mem=350GB --requeue $DIR_cur/NMTrinity.9.x.sh $nrna 

    sleep 1

done
