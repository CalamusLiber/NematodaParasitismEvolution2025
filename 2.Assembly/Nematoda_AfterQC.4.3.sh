#!/bin/sh

# module load lang/perl/5.30.0-bioperl-gcc 
# module load apps/trinity/2.8.5
module load apps/samtools/1.9
# module load apps/bowtie2/2.3.5
# module load lib/boost/1.71.0
module load lang/intel-parallel-studio-xe/2020
alias python='/To/Your/Directory/miniconda3/bin/python3.9'
alias python3='/To/Your/Directory/miniconda3/bin/python3.9'
export PATH=$PATH:/To/Your/Directory/miniconda3/bin:/To/Your/Directory/miniconda3/pkgs
# export TRINITY_HOME=/To/Your/Directory/trinityrnaseq

cd $WORK/nematoda/data/rnaseq
DIR_cur=$(pwd)
if [ ! -d trans_assembly ]; then 
    mkdir ./trans_assembly
fi

if [ ! -d AfterQC.reports ]; then 
    mkdir AfterQC.reports
fi

rnaseq_list=$(ls *_1.fastq | sed -e "s/.rna_1.fastq//g")

####Multithread preparation
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
# 将文件描述符1000与FIFO进行绑定，<读的绑定，>写的绑定，
# <>则标识对文件描述符1000的所有操作等同于对管道文件$tempfifo的操作 
rm $tmp_fifofile

for ((i=1;i<=11;i++))
do
    echo >&9
done

for nrna in $(echo $rnaseq_list | awk '{ for (i=22; i<=32; i++) {print $i} }')
do
    read -u9
    {
    fq1=$DIR_cur/$nrna.rna_1.fastq
    fq2=$DIR_cur/$nrna.rna_2.fastq

    echo "Trimming $nrna starts." >> Nema_trinity.4.log
    mkdir $nrna.AfterQC.output
    
    $HOME/pypy2/bin/pypy $HOME/AfterQC/after.py -1 $fq1 -2 $fq2 -g $nrna.AfterQC.output -b $nrna.AfterQC.output -r AfterQC.reports >> Nema_trinity.4.log 2>&1 && 
    rm $nrna.AfterQC.output/*.bad.fq && 

    echo "Trimming $nrna is done." >> Nema_trinity.4.log

    echo >&9 
    } &   
done
wait

