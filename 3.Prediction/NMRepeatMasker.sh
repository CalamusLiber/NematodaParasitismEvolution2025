#!/bin/bash -l
#SBATCH --job-name=NMsoftmask
#SBATCH --exclusive
#SBATCH --partition=skc
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --mem=350GB

cd /To/Your/Directory/nematoda/data/preds

DIR_current=$(pwd)
DIR_data=/To/Your/Directory/nematoda/data
GffRepeatMasker=$DIR_data/GffRepeatMasker.py

chmod -R 777 $GffRepeatMasker

# NNN=$(cat $DIR_current/genomic_list.dat | wc -l) # 55

tmp_fifofile="/tmp/$$.fifo" 
trap "exec 9>&-;exec 9<&-;exit 0" 2 
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$SLURM_NTASKS;i++))
do
    echo >&9
done

list=$(ls -l ./*/FilteredSeq.fna | awk -F '/' '{print $2}')
echo "Genome,FileSizeBefore,FileSizeAfter,FileSizeDiff,GenomeSizeBefore,GenomeSizeAfter,GenomeSizeDiff" > $DIR_current/summarySoftMasking.csv 

# for ((g=1;g<=$NNN;g++))
for gen in $list
do
    read -u9
    {

    # genfile=$(sed -n ${g}p $DIR_current/genomic_list.dat)
    # gen=$(echo $genfile | sed -e "s#.fna##g;s#$DIR_data\/assembly\/##g;s#_scaffolds.fasta##g;s#$DIR_data\/wgs\/genome_assembly.0\/##g")
    # genfn=$(echo $genfile | sed -e "s#$DIR_data\/assembly\/##g;s#$DIR_data\/wgs\/genome_assembly.0\/##g")
    DIR_out=$DIR_current/$gen

    mkdir -p $DIR_out 

    echo -e "#######################\nMAKER annotation start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log && 
    echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> $DIR_out/maker2pred.6.0.log && 
    time0=$(date) && start_sec=$(date --date="$time0" +%s) && 
    echo -e "#######################\nConfiguring environment...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 

    # genfile=$DIR_out/FilteredSeq.fna

    echo -e "#######################\nGenerating softmasked fasta file...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
    jump=$(grep -n '^#' $DIR_out/x01_RepeatMask.out/Repeat_elements_final_mask.gff | wc -l) && 
    $GffRepeatMasker --fasta $genfile --gff $DIR_out/x01_RepeatMask.out/Repeat_elements_final_mask.gff --out $DIR_current/$gen --skip $jump --masktype soft && 
    filesize1=$(ls -l ./$gen/FilteredSeq.fna | awk '{print $5}') 
    genomesize1=$[$(grep -v "^>" ./$gen/FilteredSeq.fna | wc -m) - $(grep -v "^>" ./$gen/FilteredSeq.fna | wc -l)]
    [[ -f ./$gen/FilteredSeq.softmask.fna ]] && filesize2=$(ls -l ./$gen/FilteredSeq.softmask.fna | awk '{print $5}') || filesize2=0 
    [[ -f ./$gen/FilteredSeq.softmask.fna ]] && genomesize2=$[$(grep -v "^>" ./$gen/FilteredSeq.softmask.fna | wc -m) - $(grep -v "^>" ./$gen/FilteredSeq.softmask.fna | wc -l)] || genomesize2=0 
     
    echo -e "$gen,$filesize1,$filesize2,$[$filesize1-$filesize2],$genomesize1,$genomesize2,$[$genomesize1-$genomesize2]" >> $DIR_current/summarySoftMasking.csv 

    echo >&9 
    } & 
done
wait

# summary

