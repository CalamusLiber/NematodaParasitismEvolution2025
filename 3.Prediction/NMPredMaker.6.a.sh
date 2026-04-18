#!/bin/bash

cd /To/Your/Directory/nematoda/data

DIR_current=$(pwd)
chmod -R 777 $DIR_current/NMPredMaker.6.a.x.sh

NNN=$(cat $DIR_current/preds/genomic_list.dat | wc -w)

for ((g=1;g<=$NNN;g++))
do

    if [[ $(ls -l $(sed -n "$g"p $DIR_current/preds/genomic_list.dat) | awk '{ print $5 }') -gt 200000000 ]]; then
        sbatch --partition=mwvdk --nodes=1 --ntasks=32 --mem=350GB --requeue $DIR_current/NMPredMaker.6.a.x.sh $g 32
    else
        sbatch --partition=short --nodes=1 --ntasks=24 --mem=100GB --requeue $DIR_current/NMPredMaker.6.a.x.sh $g 24
    fi

    sleep 1

done

# wait 

## check how many genes harvest in each run for every genome 
cd $DIR_current/preds 
echo -e "Genome,maker#1,maker#2,maker#3,augustus#1,augustus#2,snap#1,snap#2,metaeuk-ur90" > summaryMakerPred.csv 
for x in $(ls ./*/*.final.maker.pep.faa | awk -F "/" '{ print $2}') 
do 
    m[1]=$(grep '^>' ./$x/*Run1.maker.output/Run1.all.maker.proteins.fasta | wc -l)
    k=$(grep '^>' ./$x/x08_metaeuk/$x.mtek.ur90.pep.faa | wc -l)
    for i in {2..3} 
    do 
        m[$i]=$(grep '^>' ./$x/*Run$i.maker.output/Run$i.all.maker.proteins.fasta | wc -l)
        a[$i]=$(grep '^>' ./$x/*Run$i.maker.output/Run$i.all.maker.augustus.proteins.fasta | wc -l)
        s[$i]=$(grep '^>' ./$x/*Run$i.maker.output/Run$i.all.maker.snap.proteins.fasta | wc -l)
    done     
    echo "${x},${m[*]},${a[*]},${s[*]},$k" | sed -e 's/ /,/g' >> summaryMakerPred.csv 
done
