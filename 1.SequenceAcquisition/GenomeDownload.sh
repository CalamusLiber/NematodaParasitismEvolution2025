#!/bin/bash
#SBATCH --job-name=GenDown
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24Gb
#SBATCH --exclusive

# Current working directory and input CSV path.
DIR_current=$(pwd)
DIR_csv=$DIR_current/$1
SRATOOLKIT=/To/Your/Directory/sratoolkit

# Count CSV lines to determine how many rows to process.
NLine=$(cat "$DIR_csv" | wc -l)

# Resume from the last processed line if CLine.log exists.
if [ -f "$DIR_current/CLine.log" ]; then
    st=$(sed -n 1p "$DIR_current/CLine.log")
else
    # Default to line 2 to skip the CSV header.
    st=2
fi

# Create output directories for assembly, RNA-seq, and WGS data.
mkdir -p "$DIR_current/assembly"
mkdir -p "$DIR_current/rnaseq"
mkdir -p "$DIR_current/wgs"

####Multithread preparation
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=20;i++))
do
    echo >&9
done

for ((n=$st;n<=$NLine;n++))
do
    # Read the current CSV line, remove CR characters, and normalize FTP to HTTPS.
    line=$(sed -n "$n"p "$DIR_csv" | sed -e 's/\r//g;s/ftp:/https:/g')

    # Extract columns: organism, accession, type, and download URL.
    orgn=$(echo "$line" | awk -F ',' '{print $1}')
    accn=$(echo "$line" | awk -F ',' '{print $2}')
    type=$(echo "$line" | awk -F ',' '{print $3}')
    furl=$(echo "$line" | awk -F ',' '{print $4}')

    opfx="${orgn}_${accn}"

    echo "*****************"
    echo "Requesting $opfx"

    case $type in
        "Assembly_Chromosome")
            # Download chromosome-level assembly and decompress.
            wget -q -O ./assembly/$opfx.fna.gz "$furl" && \
            gzip -d ./assembly/$opfx.fna.gz
            ;;

        "Assembly_Contig" | "Assembly_Scaffold")
            # Download contig or scaffold assembly and decompress.
            wget -q -O ./assembly/$opfx.fna.gz "$furl" && \
            gzip -d ./assembly/$opfx.fna.gz
            ;;

        "RNA-Seq-TRANSCRIPTOMIC")
            # Download RNA-seq data and convert to FASTQ format.
            wget -q -O ./rnaseq/$opfx.rna "$furl"
            $SRATOOLKIT/bin/fasterq-dump -e 6 -p --split-3 --skip-technical ./rnaseq/$opfx.rna
            ;;

        "WGS-GENOMIC")
            # Download WGS genomic data and convert to FASTQ format.
            wget -q -O ./wgs/$opfx.dna "$furl"
            $SRATOOLKIT/bin/fasterq-dump -e 6 -p --split-3 --skip-technical ./wgs/$opfx.dna
            ;;
    esac

    # Wait for the current command pipeline to finish.
    wait

    # Save progress so the job can resume later if interrupted.
    echo $n > CLine.log

done
