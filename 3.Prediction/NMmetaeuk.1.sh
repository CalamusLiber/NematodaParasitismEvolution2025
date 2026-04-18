#!/bin/bash -l
#SBATCH --job-name=NMMetaEuk
#SBATCH --partition=skc
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=350GB
#SBATCH --exclusive

NP=$SLURM_CPUS_PER_TASK
# NP=$SLURM_NTASKS

DIR_cnd=/To/Your/Directory/anaconda3
DIR_current=/To/Your/Directory/nematoda/data/preds/metaeuk
DIR_data=/To/Your/Directory/nematoda/data
GffRepeatMasker=$DIR_data/GffRepeatMasker.py

chmod -R 777 $GffRepeatMasker

export PATH=/To/Your/Directory/anaconda3/bin:/To/Your/Directory/anaconda3/envs/makerenv/bin:$PATH
export LIBDIR=/To/Your/Directory/anaconda3/envs/makerenv/share/RepeatMasker/Libraries

source $DIR_cnd/bin/activate $DIR_cnd/envs/makerenv

mkdir -p $DIR_current && cd $DIR_current

cd $DIR_current

gen_list=$(ls | sed -e "s/.gen.fna//g")

tmp_fifofile="/tmp/$$.fifo" 
trap "exec 9>&-;exec 9<&-;exit 0" 2 
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$(ls *.gen.fna | wc -l);i++))
do
    echo >&9
done

for gen in $gen_list
do 
    read -u9
    {
    mkdir -p $DIR_current/$gen && 
    echo -e "#######################\nMAKER annotation start...\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log && 
    echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> $DIR_current/$gen/$gen.mtek.log && 
    time0=$(date) && start_sec=$(date --date="$time0" +%s) && 
    echo -e "#######################\nConfiguring environment...\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 

    ##---------------------Configuration---------------------##

    echo -e "#######################\nConfiguration finished.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 

    echo -e "#######################\nChecking and filtering out contigs shorter than 2 kb from genomic file.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
    if [ ! -f $DIR_current/$gen/FilteredSeq.fna ]; then 
        awk '/^>/ { if(sign_val && strLen>=2000){print sign_val ORS line} strLen=line=""; sign_val=$0; next } { strLen+=length($0); line=(line?line ORS:"")$0 } END { if(sign_val && strLen>=2000) { print sign_val ORS line }}' $DIR_current/$gen.gen.fna > $DIR_current/$gen/FilteredSeq.fna && 
        echo -e "#######################\nFiltering finished.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
        echo -e "#######################\nChecking and replacing the unsupported characters in genomic file.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
        sed -i '/^>/!y/RYKMSWBDHVrykmswbdhv/NNNNNNNNNNnnnnnnnnnn/' $DIR_current/$gen/FilteredSeq.fna && 
        echo -e "#######################\nUnsupported characters replaced by Ns.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
    else
        echo -e "#######################\nFiltering and character correction have already been done.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
    fi
    echo >&9 
    } & 
done
wait

tmp_fifofile="/tmp/$$.fifo" 
trap "exec 9>&-;exec 9<&-;exit 0" 2 
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE];i++))
do
    echo >&9
done

for gen in $gen_list
do 
    read -u9
    {
    genfile=$DIR_current/$gen/FilteredSeq.fna 
    ##---------------------Repeat Masking---------------------##
    if [ ! -f $DIR_current/$gen/Repeat_elements_final_mask.gff ]; then 
        if [ ! -f $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb/consensi.fa.classified ]; then 
            echo -e "#######################\nBuilding mask lib...\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log && 
            if [ -d $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb ]; then
                cd $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb && 
                RepeatModeler -pa $[$NP/2] -engine ncbi -database "$gen"_rpdb -dir $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb -recoverDir $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb >> $DIR_current/$gen/$gen.mtek.log 2>&1 && 
                rm -rf $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb/round-* && 
                echo -e "#######################\nMask lib built.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
            else
                mkdir -p $DIR_current/$gen/x01_RepeatMask.out && mkdir -p $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb && 
                cd $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb && 
                BuildDatabase -name "$gen"_rpdb -engine ncbi $genfile >> $DIR_current/$gen/$gen.mtek.log 2>&1 && 
                RepeatModeler -pa $[$NP/2] -engine ncbi -database "$gen"_rpdb -dir $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb >> $DIR_current/$gen/$gen.mtek.log 2>&1 && 
                rm -rf $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb/round-* && 
                echo -e "#######################\nMask lib built.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
            fi
        fi
        
        echo -e "#######################\nRepeatMasker start...\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log && 
        if [ ! -f $DIR_current/$gen/x01_RepeatMask.out/RM-specieslib/FilteredSeq.fna.masked ]; then 
            echo -e "\n#------RepeatMasking using specific lib $(date)------#\n" >> $DIR_current/$gen/$gen.mtek.log && 
            mkdir -p $DIR_current/$gen/x01_RepeatMask.out/RM-specieslib && 
            cd $DIR_current/$gen/x01_RepeatMask.out/RM-specieslib && 
            ln -s $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb/consensi.fa.classified && 
            RepeatMasker -e ncbi -pa $[$NP/2] -gff -lib $DIR_current/$gen/x01_RepeatMask.out/RM-rpdb/consensi.fa.classified -dir $DIR_current/$gen/x01_RepeatMask.out/RM-specieslib $genfile >> $DIR_current/$gen/$gen.mtek.log 2>&1 
        fi
        if [ ! -f $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib/FilteredSeq.fna.masked.out ]; then 
            echo -e "\n#------RepeatMasking using Nematoda lib from Repbase $(date)------#\n" >> $DIR_current/$gen/$gen.mtek.log && 
            mkdir -p $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib && 
            cd $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib && 
            RepeatMasker -e ncbi -pa $[$NP/2] -gff -species nematoda -dir $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib $DIR_current/$gen/x01_RepeatMask.out/RM-specieslib/FilteredSeq.fna.masked >> $DIR_current/$gen/$gen.mtek.log 2>&1 
        fi
        echo -e "\n#------Combining two libs $(date)------#\n" >> $DIR_current/$gen/$gen.mtek.log && 
        mkdir -p $DIR_current/$gen/x01_RepeatMask.out/RM-final_mask && 
        cd $DIR_current/$gen/x01_RepeatMask.out && 
        cp $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib/FilteredSeq.fna.masked.masked $DIR_current/$gen/x01_RepeatMask.out/RM-final_mask/Full_mask.fa 
        cp $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib/FilteredSeq.fna.masked.out $DIR_current/$gen/x01_RepeatMask.out/RM-final_mask/Full_mask.out 
        gunzip $DIR_current/$gen/x01_RepeatMask.out/RM-specieslib/*.cat.gz $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib/*.cat.gz 
        cat $DIR_current/$gen/x01_RepeatMask.out/RM-specieslib/*.cat $DIR_current/$gen/x01_RepeatMask.out/RM-nematodalib/*.cat > $DIR_current/$gen/x01_RepeatMask.out/RM-final_mask/Full_mask.cat 
        cd $DIR_current/$gen/x01_RepeatMask.out/RM-final_mask && 
        ProcessRepeats -species nematoda -nolow -gff $DIR_current/$gen/x01_RepeatMask.out/RM-final_mask/Full_mask.cat >> $DIR_current/$gen/$gen.mtek.log 2>&1 && 
        sed '/Satellite/d' ./Full_mask.out > Full_mask.complex.out && 
        sed '/Satellite/d' ./Full_mask.out.gff > Full_mask.complex.out.gff && 
        cat Full_mask.out.gff | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > Full_mask.reformat.gff && 
        sed 's/\t /\t/g' Full_mask.reformat.gff > $DIR_current/$gen/Repeat_elements_final_mask.gff && 
        # rm -rf $DIR_current/$gen/x01_RepeatMask.out && 
        echo -e "#######################\nRepeatMasker is done.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log
    fi
    
    if [[ ! -f $DIR_current/$gen/FilteredSeq.softmask.fna  ]]; then 
        echo -e "\n#------Generating softmasked fasta file...$(date)------#\n" >> $DIR_current/$gen/$gen.mtek.log && 
        jump=$(grep -n '^#' $DIR_current/$gen/Repeat_elements_final_mask.gff | wc -l)
        $GffRepeatMasker --fasta $genfile --gff $DIR_current/$gen/Repeat_elements_final_mask.gff --out $DIR_current/$gen --skip $jump --masktype soft 
        # HardMasker $genfile $DIR_current/$gen/Repeat_elements_final_mask.gff 
    fi

    if [ ! -f $DIR_current/$gen/$gen.mtek.ur90.pep.faa ]; then 
        echo -e "#######################\nMetaEuk (#1) prediction start...\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log && 
        mkdir -p $DIR_current/$gen/x02_metaeuk && 
        cd $DIR_current/$gen/x02_metaeuk && 
        metaeuk easy-predict ${genfile/${genfile##*.}/softmask.${genfile##*.}} /To/Your/Directory/db/UniRef90/uniref90.metaeukdb $gen.mtek.ur90 /To/Your/Directory/tmp --threads $NP --remove-tmp-files 1 --num-iterations 1 --exhaustive-search 1 --mask-lower-case 1 >> $DIR_current/$gen/$gen.mtek.log 2>&1 && 
        awk '{OFS="\t"; print $1,$3,$2,$4,$5,$6,$7,$8,$9}' $DIR_current/$gen/x02_metaeuk/$gen.mtek.ur90.gff | sed -e "s#Target_ID=#Target=#g;s#TCS_ID=#ID=#g" > $DIR_current/$gen/$gen.mtek.ur90.reform.gff && 
        mv $gen.mtek.ur90.fas $DIR_current/$gen/$gen.mtek.ur90.pep.faa && cp -f $DIR_current/$gen/$gen.mtek.ur90.pep.faa $DIR_data/busco/pep && 
        mv $gen.mtek.ur90.codon.fas $DIR_current/$gen/$gen.mtek.ur90.cds.fna && # cp -f $gen.mtek.ur90.cds.fna $DIR_data/busco/cds && 
        rm -rf $DIR_current/$gen/x02_metaeuk && 
        echo -e "#######################\nMetaEuk (#1) prediction is done.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
    else
        echo -e "#######################\nThe file $DIR_current/$gen/$gen.mtek.ur90.pep.faa is already there.\n$(date)\n#######################\n\n" >> $DIR_current/$gen/$gen.mtek.log 
    fi
    echo >&9 
    } & 
done
wait
