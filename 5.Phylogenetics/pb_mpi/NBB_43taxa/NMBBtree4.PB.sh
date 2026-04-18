#!/bin/bash

module load apps/bmge/1.12/

cd /To/Your/Directory/nematoda/phylo/AllBackbone4.nom.rmsp

BMGE=/sw/apps/BMGE-1.12/BMGE.jar

echo -e "=======================\nTrimming with BMGE (normal) start...\n$(date)\n=======================\n\n" >> BBmsa.bmge.log && 
java -Xmx90G -jar $BMGE -i ./NBBData43taxa174426AA.phy -t AA -m BLOSUM90 -h 0.4 -s NO -op ./NBBData43taxa174426AA.bmge.phy >> BBmsa.bmge.log 2>&1 && 
echo -e "=======================\nTrimming with BMGE (normal) accomplished.\n$(date)\n=======================\n\n" >> BBmsa.bmge.log

chmod -R 777 ./NMPB4.x.sh

for run in {1..6} 
do 
    sbatch --partition=compute ./NMPB4.x.sh $run 
done

## for the second, third, or more 2000 generations, please use the following code:
# for run in 1 3 5 # change run number according to the actual situation 
# do 
#     sed -i 's#1\t6000#1\t8000#g' AllBackbone4.nom.rmsp.pb.chain$run.param  # change generations according to the actual situation
# 	  sbatch ./NMPB4.x.sh $run 
# done


