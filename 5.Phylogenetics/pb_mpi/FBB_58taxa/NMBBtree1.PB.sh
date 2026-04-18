#!/bin/bash

DIR_cnd=/public1/home/scb4616/apps/anaconda3
source $DIR_cnd/bin/activate phylobayes

cd /public1/home/scb4616/AllBackbone1.nom.rmsp

BMGE=/sw/apps/BMGE-1.12/BMGE.jar

echo -e "=======================\nTrimming with BMGE (normal) start...\n$(date)\n=======================\n\n" >> BBmsa.bmge.log && 
java -Xmx90G -jar $BMGE -i ./FBBData58taxa149581AA.phy -t AA -m BLOSUM45 -h 0.3 -s NO -op ./FBBData58taxa149581AA.bmge.phy >> BBmsa.bmge.log 2>&1 && 
echo -e "=======================\nTrimming with BMGE (normal) accomplished.\n$(date)\n=======================\n\n" >> BBmsa.bmge.log

# chmod -R 777 ./NMPB1.x.sh

for run in {1..6} 
do 
    sbatch ./NMPB1.x.sh $run 
done

## for the second, third, or more 2000 generations, please use the following code:
# for run in 1 3 5 # change run number according to the actual situation 
# do 
#     sed -i 's#1\t6000#1\t8000#g' AllBackbone1.nom.rmsp.pb.chain$run.param  # change generations according to the actual situation
# 	sbatch ./NMPB1.x.sh $run 
# done


