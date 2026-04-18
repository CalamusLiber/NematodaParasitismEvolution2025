#!/bin/sh
#SBATCH --job-name=NMTrinity.9

cd /To/Your/Directory/nematoda/data/rnaseq
DIR_cur=$(pwd)
DIR_cnd=/To/Your/Directory/anaconda3
DIR_asm=$DIR_cur/trans_assembly

source $DIR_cnd/bin/activate $DIR_cnd/envs/rnaseqenv

rnaseq_list=$(ls $DIR_asm/*.Trinity.fasta | sed -e "s#$DIR_asm\/##g;s#.Trinity.fasta##g")

echo -e "#######################\nTrinity Assembling Quality Assessment...\n$(date)\n#######################\n\n" >> $DIR_asm/NMTrinityQA.log && 
echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\nNumber of CPUS: $SLURM_NTASKS\n" >> $DIR_asm/NMTrinityQA.log && 

tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile
for ((i=1;i<=$SLURM_NTASKS;i++))
do
    echo >&9
done

for nrna in $rnaseq_list 
do 
    read -u9
    {

    if [ ! -f $DIR_asm/TrinityAssemblingStat.csv ]; then
        echo "Assembly,TotalGenes,TotalTranscripts,%GC,All-ContigN50,All-MedianLength,All-AveLength,All-TotalBases,LLL-ContigN50,LLL-MedianLength,LLL-AveLength,LLL-TotalBases" > $DIR_asm/TrinityAssemblingStat.csv  
    fi

    qares=$($DIR_cnd/envs/rnaseqenv/bin/TrinityStats.pl $DIR_asm/$nrna.Trinity.fasta) && 
    qaline=$(echo $qares | sed -e "s/################################ ## Counts of transcripts, etc. ################################ Total trinity 'genes'://g;s/Total trinity transcripts://g;s/Percent GC://g;s/######################################## Stats based on ALL transcript contigs: ######################################## //g;s/Contig N10://g;s/Contig N20://g;s/Contig N30://g;s/Contig N40://g;s/Contig N50://g;s/Median contig length://g;s/Average contig://g;s/Total assembled bases://g;s/##################################################### ## Stats based on ONLY LONGEST ISOFORM per 'GENE': #####################################################//g") && 
    
    echo $nrna,$(echo $qaline | awk '{ print $1, $2, $3, $8, $9, $10, $11, $16, $17, $18, $19 }' | sed "s/ /,/g") >> $DIR_asm/TrinityAssemblingStat.csv 
    
    echo >&9 
    } & 
done

wait

echo -e "#######################\nTrinity Assembling Quality Assessment is finished.\n#######################\n\n" >> $DIR_asm/NMTrinityQA.log 
