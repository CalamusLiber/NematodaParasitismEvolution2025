#!/bin/bash

DIR_cur=$(pwd)
cd $DIR_cur

Topo=$(echo $(pwd) | awk -v FS='/' '{print $NF}' | awk -v FS='.BD' '{print $1}')
TreeCalibration=$(ls *.tre | sed "s#\/##g;s#Nema.${Topo}.##g;s#.tre##g")
ClockModel='IR AR'
Partition='1P 5P'
rgene_gamma='15 20 25'

# Run MCMCtree step 1 script
for TreCali in $TreeCalibration
do 
    for Part in $Partition 
    do
        for Clock in $ClockModel
        do 
            for Rgene in $rgene_gamma
            do 
                if [ ! -f $DIR_cur/${TreCali}.${Part}.${Clock}.rg${Rgene}/*MCMCtree_final_report.txt ] 
                then 
                    sbatch $DIR_cur/NMMCMCTr.8.x.1.sh $TreCali $Part $Clock $Rgene
                fi
            done
        done
    done
done

