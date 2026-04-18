#!/bin/sh

cd /To/Your/Directory/nematoda/data

for k in 15 19 23 27 31 35 39 43 47 
do
    sbatch --job-name=NMSPT$k --partition=short --nodes=1 --ntasks=24 --mem=150Gb --exclusive ./NMSOAP_trans.x.sh $k 
done

# run the following commands after all the NMSOAP_trans.x.sh being finished
findK(){
    local n=1 
    local Kmer=(15 19 23 27 31 35 39 43 47)
    for i in $1 
    do 
        if [[ $i -eq $2 ]] 
        then
            echo Best-K=${Kmer[$n-1]}
        fi
        ((n++))
    done
}

Kmer=(15 19 23 27 31 35 39 43 47)

scaffoldN50=$(grep -P 'N50\t' ./*scafStatistics | sed -n '1~2p' | awk '{print $2}')
scaffoldN50Max=$(grep -P 'N50\t' ./*scafStatistics | sed -n '1~2p' | awk '{print $2}' | sort -rn -t ' ' | head -1)
scaffoldMedSize=$(grep -P 'Median_Size\t' ./*scafStatistics | sed -n '1~2p' | awk '{print $2}')
scaffoldMedSizeMax=$(grep -P 'Median_Size\t' ./*scafStatistics | sed -n '1~2p' | awk '{print $2}' | sort -rn -t ' ' | head -1)
scaffoldMeanSize=$(grep -P 'Mean_Size\t' ./*scafStatistics | sed -n '1~2p' | awk '{print $2}')
scaffoldMeanSizeMax=$(grep -P 'Mean_Size\t' ./*scafStatistics | sed -n '1~2p' | awk '{print $2}' | sort -rn -t ' ' | head -1)
scaffold1k=$(grep 'scaffolds>1K' ./*scafStatistics | awk '{print $2}')
scaffold1kMax=$(grep 'scaffolds>1K' ./*scafStatistics | awk '{print $2}' | sort -rn -t ' ' | head -1)

findK "$scaffoldN50" $scaffoldN50Max # Best-K=23
findK "$scaffoldMedSize" $scaffoldMedSizeMax # Best-K=43
findK "$scaffoldMeanSize" $scaffoldMeanSizeMax # Best-K=43
findK "$scaffold1k" $scaffold1kMax # Best-K=23

contigN50=$(grep -P 'N50\t' ./*scafStatistics | sed -n '0~2p' | awk '{print $2}')
contigN50Max=$(grep -P 'N50\t' ./*scafStatistics | sed -n '0~2p' | awk '{print $2}' | sort -rn -t ' ' | head -1)
contigMedSize=$(grep -P 'Median_Size\t' ./*scafStatistics | sed -n '0~2p' | awk '{print $2}')
contigMedSizeMax=$(grep -P 'Median_Size\t' ./*scafStatistics | sed -n '0~2p' | awk '{print $2}' | sort -rn -t ' ' | head -1)
contigMeanSize=$(grep -P 'Mean_Size\t' ./*scafStatistics | sed -n '0~2p' | awk '{print $2}')
contigMeanSizeMax=$(grep -P 'Mean_Size\t' ./*scafStatistics | sed -n '0~2p' | awk '{print $2}' | sort -rn -t ' ' | head -1)
contig500=$(grep 'Contig>500' ./*scafStatistics | awk '{print $2}')
contig500Max=$(grep 'Contig>500' ./*scafStatistics | awk '{print $2}' | sort -rn -t ' ' | head -1)

findK "$contigN50" $contigN50Max # Best-K=19 
findK "$contigMedSize" $contigMedSizeMax # Best-K=47 
findK "$contigMeanSize" $contigMeanSizeMax # Best-K=19 
findK "$contig500" $contig500Max # Best-K=19 

# select scaffolf of K=23 as the best assembly by SOAPdenovo-trans (suffix with _trinity.Trinity is used only for handling by subsequent commmands)
cp /To/Your/Directory/nematoda/data/rnaseq/Enoplus_brevis_NEB.A64176.soap_trans/res.k23.scafSeq /To/Your/Directory/nematoda/data/rnaseq/trans_assembly/Enoplus_brevis_NEB.A64176.soap.fasta 

