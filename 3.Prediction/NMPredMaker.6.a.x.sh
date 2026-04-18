#!/bin/bash -l
#SBATCH --job-name=NMPM.x
#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=48
#SBATCH --mem=360GB
#SBATCH --exclusive

NP=$2

module load lib/openmpi/4.0.2-gcc.4.8.5

cd /To/Your/Directory/nematoda/data/preds

DIR_current=$(pwd)
DIR_data=/To/Your/Directory/nematoda/data
DIR_cnd=/To/Your/Directory/anaconda3
GffRepeatMasker=$DIR_data/GffRepeatMasker.py

chmod -R 777 $GffRepeatMasker

genfile=$(sed -n "$1"p $DIR_current/genomic_list.dat)
gen=$(echo $genfile | sed -e "s#.fna##g;s#$DIR_data\/assembly\/##g;s#_scaffolds.fasta##g;s#$DIR_data\/wgs\/genome_assembly.0\/##g")
genfn=$(echo $genfile | sed -e "s#$DIR_data\/assembly\/##g;s#$DIR_data\/wgs\/genome_assembly.0\/##g")
DIR_out=$DIR_current/$gen

mkdir -p $DIR_out 

echo -e "#######################\nMAKER annotation start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log && 
echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> $DIR_out/maker2pred.6.0.log && 
time0=$(date) && start_sec=$(date --date="$time0" +%s) && 
echo -e "#######################\nConfiguring environment...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 

export PATH=/To/Your/Directory/anaconda3/bin:/To/Your/Directory/anaconda3/envs/makerenv/bin:$PATH
export LIBDIR=/To/Your/Directory/anaconda3/envs/makerenv/share/RepeatMasker/Libraries
# export LD_PRELOAD=/To/Your/Directory/anaconda3/pkgs/mpich-3.4.2-h846660c_100/lib/libmpi.so
# export LD_PRELOAD=/sw/lib/openmpi-4.0.2-gcc-4.8.5/lib/libmpi.so

source $DIR_cnd/bin/activate $DIR_cnd/envs/makerenv

est_evid=$(echo $(ls $DIR_current/ests/*.fa)" "$(ls $DIR_current/ests/*.fna)" "$(ls $DIR_current/ests/*.fasta) | sed -e 's/ /,/g')
prot_evid=$(echo $(ls $DIR_current/proteins/*.faa)" "$(ls $DIR_current/proteins/*.protein.fa) | sed -e 's/ /,/g')",$DIR_current/proteins/uniprot_Ecdysozoa_reviewed_canonical+isoform.fasta,$DIR_current/proteins/GeneDatabaseNematoda20211201112533.pept.fasta"

##---------------------Function---------------------##
trainAugustus(){
    if [ ! -d $1 ]; then
        mkdir -p $1 && mkdir -p $1/species && 
        cp -r $2/cgp $1 && 
        cp -r $2/extrinsic $1 && 
        cp -r $2/model $1 && 
        cp -r $2/profile $1 && 
        cp -r $2/species/generic $1/species && 
        cp -r $2/species/ancylostoma_ceylanicum $1/species && 
        cp -r $2/species/brugia $1/species && 
        cp -r $2/species/caenorhabditis $1/species && 
        cp -r $2/species/c_elegans_trsk $1/species && 
        cp -r $2/species/trichinella $1/species 
    fi
    export AUGUSTUS_CONFIG_PATH=$1 && 
    awk '{if ($2=="maker") print }' $5 > $8.maker.gff && 
    gff2gbSmallDNA.pl $8.maker.gff $6 1000 $8.maker.gb && ## We keep 1000 bp up- and down-stream of each gene for training the models based on our initial MAKER run
    echo -e "-----------------------\n$(grep -c LOCUS $8.maker.gb) genes in training set.\n-----------------------\n\n" && 
    echo -e "-----------------------\nSplitting gene structure set into a training and a test set\n$(date)\n-----------------------\n\n" && 
    testnum=$(echo "scale=0; $(grep -c LOCUS $8.maker.gb) * $3 / 1" | bc) && 
    echo -e "Test set has $testnum genes.\n" && 
    randomSplit.pl $8.maker.gb $testnum && 
    grep -c LOCUS $8.maker.gb* && 
    echo -e "\n-----------------------\n\n" && 
    echo -e "-----------------------\nCreate a meta parameters file for current species\n$(date)\n-----------------------\n\n" && 
    new_species.pl --species=$7.$8 && 
    echo -e "-----------------------\nInitial training\n$(date)\n-----------------------\n\n" && 
    etraining --species=$7.$8 $8.maker.gb && 
    ls -ort $AUGUSTUS_CONFIG_PATH/species/$7.$8 && 
    echo -e "-----------------------\nEvaluation before optimization\n$(date)\n" && 
    augustus --species=$7.$8 $8.maker.gb.test >& $8.maker.gb.test1.out && 
    grep -A 22 Evaluation $8.maker.gb.test1.out && 
    echo -e "\n-----------------------\n\n" && 
    echo -e "-----------------------\nOptimizes the prediction accuracy by adjusting the meta parameters by $4 rounds\n$(date)\n-----------------------\n\n" && 
    optimize_augustus.pl --species=$7.$8 --kfold=$9 --cpus=$9 --rounds=$4 --cleanup=1 --onlytrain=$8.maker.gb.train $8.maker.gb.test && 
    echo -e "-----------------------\nPost-optimization training\n$(date)\n-----------------------\n\n" && 
    etraining --species=$7.$8 $8.maker.gb && 
    echo -e "-----------------------\nEvaluation after optimization\n$(date)\n" && 
    augustus --species=$7.$8 $8.maker.gb.test >& $8.maker.gb.test2.out && 
    grep -A 22 Evaluation $8.maker.gb.test2.out 
}

HardMasker(){
    local infna=$1 
    local gff=$2 
    local outfna=${infna/${infna##*.}/hardmask.${infna##*.}} 
    
    # get fna line numbers and contig tags 
    local LineNum=($(grep -n "^>" $infna | awk -F ":>" '{print $1}')) 
    local ContigTagQ=($(grep "^>" $infna | awk -F "[> ]" '{print $2}')) 
    
    # convert sequence into a contig array on which we change the mask positions into N 
    for g in ${!ContigTagQ[*]} 
    do 
        # get the seq of of the #g th contig and convert to an array
        local sline=$[${LineNum[$g]}+1] 
        local eline=$[${LineNum[$[$g+1]]}-1] 
        local seqx=($(sed -n "${sline},${eline}p" $infna | sed ':t;N;s/\n//;b t' | sed -e "s/./& /g")) 
        
        # get substitute positions of the g th contig
        local Qstart=($(grep "^${ContigTagQ[$g]}" $gff | awk '{ print $4 }')) 
        local Qend=($(grep "^${ContigTagQ[$g]}" $gff | awk '{ print $5 }')) 

        local xpos=()
        for i in ${!Qstart[*]}
        do 
            xpos=(${xpos[*]} $(seq $[${Qstart[$i]}-1] $[${Qend[$i]}-1])) 
        done
        
        # replace the bases at xpos with N (parallel execution)
        for j in ${xpos[*]} 
        do 
            seqx[$j]=N 
        done

        seqx=$(echo "${seqx[*]}" | sed "s/ //g") 

        # write the masked seq into local disc (80 bases per line)
        sed -n "${LineNum[$g]}p" $infna >> $outfna 
        local len=${#seqx} 
        local nline=$((len/80)) 
        local ntail=$((len%80)) 
        local sp=0 
        for ((k=1;k<=$nline;k++)) 
        do 
            echo ${seqx:sp:80} >> $outfna 
            ((sp+=80))
        done
        [ $ntail -ne 0 ] && echo ${seqx:sp:ntail} >> $outfna 
    done
}

##---------------------Configuration---------------------##

echo -e "#######################\nConfiguration finished.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 

echo -e "#######################\nChecking and filtering out contigs shorter than 2 kb from genomic file.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 
if [ ! -f $DIR_out/FilteredSeq.fna ]; then 
    awk '/^>/ { if(sign_val && strLen>=2000){print sign_val ORS line} strLen=line=""; sign_val=$0; next } { strLen+=length($0); line=(line?line ORS:"")$0 } END { if(sign_val && strLen>=2000) { print sign_val ORS line }}' $genfile > $DIR_out/FilteredSeq.fna && 
    genfile=$DIR_out/FilteredSeq.fna && 
    echo -e "#######################\nFiltering finished.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 
    echo -e "#######################\nChecking and replacing the unsupported characters in genomic file.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 
    sed -i '/^>/!y/RYKMSWBDHVrykmswbdhv/NNNNNNNNNNnnnnnnnnnn/' $genfile && 
    echo -e "#######################\nUnsupported characters replaced by Ns.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 
else
    genfile=$DIR_out/FilteredSeq.fna && 
    echo -e "#######################\nFiltering and character correction have already been done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 
fi

##---------------------Repeat Masking---------------------##
if [ ! -f $DIR_out/x01_RepeatMask.out/Repeat_elements_final_mask.gff ]; then 
    if [ ! -f $DIR_out/x01_RepeatMask.out/RM-rpdb/consensi.fa.classified ]; then 
        echo -e "#######################\nBuilding mask lib...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log && 
        if [ -d $DIR_out/x01_RepeatMask.out/RM-rpdb ]; then
            cd $DIR_out/x01_RepeatMask.out/RM-rpdb && 
            RepeatModeler -pa $[$NP-1] -engine ncbi -database "$gen"_rpdb -dir $DIR_out/x01_RepeatMask.out/RM-rpdb -recoverDir $DIR_out/x01_RepeatMask.out/RM-rpdb >> $DIR_out/maker2pred.6.0.log 2>&1 && 
            rm -rf $DIR_out/x01_RepeatMask.out/RM-rpdb/round-* && 
            echo -e "#######################\nMask lib built.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 
        else
            mkdir -p $DIR_out/x01_RepeatMask.out && mkdir -p $DIR_out/x01_RepeatMask.out/RM-rpdb && 
            cd $DIR_out/x01_RepeatMask.out/RM-rpdb && 
            BuildDatabase -name "$gen"_rpdb -engine ncbi $genfile >> $DIR_out/maker2pred.6.0.log 2>&1 && 
            RepeatModeler -pa $[$NP-1] -engine ncbi -database "$gen"_rpdb -dir $DIR_out/x01_RepeatMask.out/RM-rpdb >> $DIR_out/maker2pred.6.0.log 2>&1 && 
            rm -rf $DIR_out/x01_RepeatMask.out/RM-rpdb/round-* && 
            echo -e "#######################\nMask lib built.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log 
        fi
    fi
    
    echo -e "#######################\nRepeatMasker start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log && 
    if [ ! -f $DIR_out/x01_RepeatMask.out/RM-specieslib/FilteredSeq.fna.masked ]; then 
        echo -e "\n#------RepeatMasking using specific lib $(date)------#\n" >> $DIR_out/maker2pred.6.0.log && 
        mkdir -p $DIR_out/x01_RepeatMask.out/RM-specieslib && 
        cd $DIR_out/x01_RepeatMask.out/RM-specieslib && 
        ln -s $DIR_out/x01_RepeatMask.out/RM-rpdb/consensi.fa.classified && 
        RepeatMasker -e ncbi -pa $[$NP-1] -gff -lib $DIR_out/x01_RepeatMask.out/RM-rpdb/consensi.fa.classified -dir $DIR_out/x01_RepeatMask.out/RM-specieslib $genfile >> $DIR_out/maker2pred.6.0.log 2>&1 
    fi
    if [ ! -f $DIR_out/x01_RepeatMask.out/RM-nematodalib/FilteredSeq.fna.masked.out ]; then 
        echo -e "\n#------RepeatMasking using Nematoda lib from Repbase $(date)------#\n" >> $DIR_out/maker2pred.6.0.log && 
        mkdir -p $DIR_out/x01_RepeatMask.out/RM-nematodalib && 
        cd $DIR_out/x01_RepeatMask.out/RM-nematodalib && 
        RepeatMasker -e ncbi -pa $[$NP-1] -gff -species nematoda -dir $DIR_out/x01_RepeatMask.out/RM-nematodalib $DIR_out/x01_RepeatMask.out/RM-specieslib/FilteredSeq.fna.masked >> $DIR_out/maker2pred.6.0.log 2>&1 
    fi
    echo -e "\n#------Combining two libs $(date)------#\n" >> $DIR_out/maker2pred.6.0.log && 
    mkdir -p $DIR_out/x01_RepeatMask.out/RM-final_mask && 
    cd $DIR_out/x01_RepeatMask.out && 
    cp $DIR_out/x01_RepeatMask.out/RM-nematodalib/FilteredSeq.fna.masked.masked $DIR_out/x01_RepeatMask.out/RM-final_mask/Full_mask.fa 
    cp $DIR_out/x01_RepeatMask.out/RM-nematodalib/FilteredSeq.fna.masked.out $DIR_out/x01_RepeatMask.out/RM-final_mask/Full_mask.out 
    gunzip $DIR_out/x01_RepeatMask.out/RM-specieslib/*.cat.gz $DIR_out/x01_RepeatMask.out/RM-nematodalib/*.cat.gz 
    cat $DIR_out/x01_RepeatMask.out/RM-specieslib/*.cat $DIR_out/x01_RepeatMask.out/RM-nematodalib/*.cat > $DIR_out/x01_RepeatMask.out/RM-final_mask/Full_mask.cat 
    cd RM-final_mask && 
    ProcessRepeats -species nematoda -nolow -gff $DIR_out/x01_RepeatMask.out/RM-final_mask/Full_mask.cat >> $DIR_out/maker2pred.6.0.log 2>&1 && 
    sed '/Satellite/d' ./Full_mask.out > Full_mask.complex.out && 
    sed '/Satellite/d' ./Full_mask.out.gff > Full_mask.complex.out.gff && 
    cat Full_mask.out.gff | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > Full_mask.reformat.gff && 
    sed 's/\t /\t/g' Full_mask.reformat.gff > $DIR_out/x01_RepeatMask.out/Repeat_elements_final_mask.gff && 
    echo -e "#######################\nRepeatMasker is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.0.log
fi

##---------------------MAKER RUN #1---------------------##
# Use "sed -i '2~2y/RYKMSWBDHV/NNNNNNNNNN/' GeneDatabaseNematoda20211201112533.mRNA.fasta" to replace unsupported characters with N or add -fix_necleotides in the command when they were checked out and block the execution of maker.
if [ ! -f $DIR_out/x02_Run1.maker.output/Run1.all.gff ]; then
    if [ -f $DIR_out/x01_RepeatMask.out/Repeat_elements_final_mask.gff ]; then # $DIR_out/x01_RepeatMask.out/RM-rpdb/consensi.fa.classified

        cp $DIR_current/maker*.ctl $DIR_out && 
        cd $DIR_out && 

        echo -e "#######################\nChange configuration file for MAKER (run #1) start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log && 
        sed -i -r "s#^genome=#&$genfile#g" $DIR_out/maker_opts.ctl && 
        # sed -i -e "s#^est=#&$DIR_current\/ests\/GeneDatabaseNematoda20211201112533.mRNA.fasta#g" $DIR_out/maker_opts.ctl && 
        sed -i -e "s#^est=#&$est_evid#g" $DIR_out/maker_opts.ctl && 
        # sed -i -e "s#^protein=#&$DIR_current\/proteins\/uniprot_Ecdysozoa_reviewed_canonical+isoform.fasta,$DIR_current\/proteins\/GeneDatabaseNematoda20211201112533.pept.fasta#g" $DIR_out/maker_opts.ctl && 
        sed -i -e "s#^protein=#&$prot_evid#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^model_org=#&simple#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^rmlib=#&$DIR_out\/x01_RepeatMask.out\/RM-rpdb\/consensi.fa.classified#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^rm_gff=#&$DIR_out\/x01_RepeatMask.out\/Repeat_elements_final_mask.gff#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^alt_splice=0#alt_splice=1#g" $DIR_out/maker_opts.ctl && 
        echo -e "#######################\nChange configuration file for MAKER (run #1) processes are done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log 

        echo -e "#######################\nMAKER (run #1) start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log && 
        /sw/lib/openmpi-4.0.2-gcc-4.8.5/bin/mpiexec -np $NP --oversubscribe maker --ignore_nfs_tmp -base Run1 >> $DIR_out/maker2pred.6.1.log 2>&1 && 
        mv $DIR_out/Run1.maker.output $DIR_out/x02_Run1.maker.output && 
        echo -e "#######################\nMAKER (run #1) processes are done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log 

        if [[ $(ls -l $DIR_out/x02_Run1.maker.output/Run1_master_datastore_index.log | awk '{ print $5 }') -gt 0 ]]; then

            cd $DIR_out/x02_Run1.maker.output && 

            echo -e "#######################\nCombining MAKER (run #1) results...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log && 

            gff3_merge -d $DIR_out/x02_Run1.maker.output/Run1_master_datastore_index.log && 
            fasta_merge -d $DIR_out/x02_Run1.maker.output/Run1_master_datastore_index.log && 
            rm -rf $DIR_out/x02_Run1.maker.output/mpi_blastdb && 
            rm -rf $DIR_out/x02_Run1.maker.output/Run1_datastore && 
            rm -f $DIR_out/x02_Run1.maker.output/.NFSLock.* 
            rm -rf $DIR_out/x02_Run1.maker.output/repeatmasker.first.* 
            
            echo -e "#######################\nMAKER (run #1) job is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log 
        else
            echo -e "#######################\nMAKER (run #1) job is finished but with no gain.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log  
        fi
    else
        echo -e "#######################\nNo RepeatMasking results found.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.1.log  
    fi
fi

# source $DIR_cnd/bin/deactivate

##---------------------MAKER RUN #2---------------------##
if [ ! -f $DIR_out/x05_Run2.maker.output/Run2.all.gff ]; then
    if [ -f $DIR_out/x02_Run1.maker.output/Run1.all.gff ]; then
        
        if [ ! -f $DIR_out/x03_snap1/$gen.snap1.hmm ]; then
            echo -e "#######################\nSNAP (#1) training start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log && 
            cd $DIR_out/x02_Run1.maker.output && 
            maker2zff -n Run1.all.gff && 
            mkdir -p $DIR_out/x03_snap1 && 
            cp $DIR_out/x02_Run1.maker.output/genome.* $DIR_out/x03_snap1 && 
            cd $DIR_out/x03_snap1 && 
            fathom -categorize 1000 genome.ann genome.dna && 
            fathom -export 1000 -plus uni.ann uni.dna && 
            forge export.ann export.dna && 
            hmm-assembler.pl $gen.snap1 . > $gen.snap1.hmm && 
            echo -e "#######################\nSNAP (#1) training is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log  
        else
            echo -e "#######################\nThe file $DIR_out/x03_snap1/$gen.snap1.hmm is already there.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log 
        fi
        
        if [ ! -d $DIR_current/augustus_config/species/$gen.Run1 ]; then
            echo -e "#######################\nAugustus (#1) training start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log && 
            mkdir -p $DIR_out/x04_augustus1 && 
            cd $DIR_out/x04_augustus1 && 
            trainAugustus $DIR_current/augustus_config $DIR_cnd/envs/makerenv/config 0.25 5 $DIR_out/x02_Run1.maker.output/Run1.all.gff $genfile $gen Run1 $NP >> $DIR_out/maker2pred.6.2.log 2>&1 && 
            cp -r $DIR_current/augustus_config/species/$gen.Run1 . && 
            echo -e "#######################\nAugustus (#1) training is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log  
        else
            export AUGUSTUS_CONFIG_PATH=$DIR_current/augustus_config && 
            echo -e "#######################\nThe model $gen.Run1 is ready.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log 
        fi
        
        echo -e "#######################\nGeneMark-ES (#1) training start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log && 
        echo -e "Currently not provided, as we do not have the mRANA-Seq data that have been generated as part of the original sequencing project." >> $DIR_out/maker2pred.6.2.log && 
        echo -e "#######################\nGeneMark-ES (#1) training is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log && 

        cp $DIR_current/maker*.ctl $DIR_out &&
        cd $DIR_out && 

        echo -e "#######################\nChange configuration file for MAKER (run #2) start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log && 
        sed -i -r "s#^genome=#&$genfile#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^maker_gff=#&$DIR_out\/x02_Run1.maker.output\/Run1.all.gff#g" $DIR_out/maker_opts.ctl && 
        # sed -i -e "s#^protein=#&$DIR_current\/proteins\/uniprot_Ecdysozoa_reviewed_canonical+isoform.fasta#g;s#^protein_gff=#&$DIR_current\/proteins\/uniprot_Ecdysozoa_reviewed.gff#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^est_pass=0#est_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^altest_pass=0#altest_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^protein_pass=0#protein_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^rm_pass=0#rm_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^model_pass=0#model_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^pred_pass=0#pred_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^other_pass=0#other_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^est2genome=1#est2genome=0#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^protein2genome=1#protein2genome=0#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^snaphmm=#&$DIR_out\/x03_snap1\/$gen.snap1.hmm#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^gmhmm=#&$DIR_out\/genemarkes1\/$gen.gm1.hmm#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^augustus_species=#&$gen.Run1#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^run_evm=0#run_evm=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^pred_stats=0#pred_stats=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^alt_splice=0#alt_splice=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^keep_preds=0#keep_preds=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^correct_est_fusion=1#correct_est_fusion=0#g" $DIR_out/maker_opts.ctl && 
        echo -e "#######################\nChange configuration file for MAKER (run #2) processes are done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log 

        echo -e "#######################\nMAKER (run #2) start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log && 
        /sw/lib/openmpi-4.0.2-gcc-4.8.5/bin/mpiexec -np $NP --oversubscribe maker --ignore_nfs_tmp -base Run2 -RM_off >> $DIR_out/maker2pred.6.2.log 2>&1 && 
        mv $DIR_out/Run2.maker.output $DIR_out/x05_Run2.maker.output && 
        echo -e "#######################\nMAKER (run #2) processes are done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log 

        if [[ $(ls -l $DIR_out/x05_Run2.maker.output/Run2_master_datastore_index.log | awk '{ print $5 }') -gt 0 ]]; then

            cd $DIR_out/x05_Run2.maker.output && 

            echo -e "#######################\nCombining MAKER (run #2) results...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log && 

            gff3_merge -d $DIR_out/x05_Run2.maker.output/Run2_master_datastore_index.log && 
            fasta_merge -d $DIR_out/x05_Run2.maker.output/Run2_master_datastore_index.log && 
            rm -rf $DIR_out/x05_Run2.maker.output/mpi_blastdb && 
            rm -rf $DIR_out/x05_Run2.maker.output/Run2_datastore && 
            rm -f $DIR_out/x05_Run2.maker.output/.NFSLock.* 

            echo -e "#######################\nMAKER (run #2) job is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log 
        else
            echo -e "#######################\nMAKER (run #2) job is finished but with no gain.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log  
        fi
    else
        echo -e "#######################\nNo Run1.all.gff found...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.2.log 
    fi
fi

# source $DIR_cnd/bin/deactivate

##---------------------MAKER RUN #3---------------------##
if [ ! -f $DIR_out/x09_Run3.maker.output/Run3.all.gff ]; then
    if [ -f $DIR_out/x05_Run2.maker.output/Run2.all.gff ]; then
        
        if [ ! -f $DIR_out/x06_snap2/$gen.snap2.hmm ]; then
            echo -e "#######################\nSNAP (#2) training start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
            cd $DIR_out/x05_Run2.maker.output && 
            maker2zff -n Run2.all.gff && 
            mkdir -p $DIR_out/x06_snap2 && 
            mv $DIR_out/x05_Run2.maker.output/genome.* $DIR_out/x06_snap2 && 
            cd $DIR_out/x06_snap2 && 
            fathom -categorize 1000 genome.ann genome.dna && 
            fathom -export 1000 -plus uni.ann uni.dna && 
            forge export.ann export.dna && 
            hmm-assembler.pl $gen.snap2 . > $gen.snap2.hmm && 
            echo -e "#######################\nSNAP (#2) training is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
        else
            echo -e "#######################\nThe file $DIR_out/x06_snap2/$gen.snap2.hmm is already there.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
        fi
        
        if [ ! -d $DIR_current/augustus_config/species/$gen.Run2 ]; then
            echo -e "#######################\nAugustus (#2) training start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
            mkdir -p $DIR_out/x07_augustus2 && 
            cd $DIR_out/x07_augustus2 && 
            trainAugustus $DIR_current/augustus_config $DIR_cnd/envs/makerenv/config 0.25 5 $DIR_out/x05_Run2.maker.output/Run2.all.gff $genfile $gen Run2 $NP >> $DIR_out/maker2pred.6.3.log 2>&1 && 
            cp -r $DIR_current/augustus_config/species/$gen.Run2 . && 
            echo -e "#######################\nAugustus (#2) training is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
        else
            export AUGUSTUS_CONFIG_PATH=$DIR_current/augustus_config && 
            echo -e "#######################\nThe model $gen.Run2 is ready.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
        fi

        echo -e "#######################\nGeneMark-ES (#2) training start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
        echo -e "Currently not provided, as we do not have the mRNA-Seq data that have been generated as part of the original sequencing project." >> $DIR_out/maker2pred.6.3.log && 
        echo -e "#######################\nGeneMark-ES (#2) training is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 

        if [ ! -f $DIR_out/x08_metaeuk/$gen.mtek.ur90.pep.faa ]; then 
            echo -e "#######################\nGenerating hardmasked fasta file...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
            jump=$(grep -n '^#' $DIR_out/x01_RepeatMask.out/Repeat_elements_final_mask.gff | wc -l) && 
            $GffRepeatMasker --fasta $genfile --gff $DIR_out/x01_RepeatMask.out/Repeat_elements_final_mask.gff --out $DIR_current/$gen --skip $jump --masktype soft && 
            echo -e "#######################\nMetaEuk (#1) prediction start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
            mkdir -p $DIR_out/x08_metaeuk && 
            cd $DIR_out/x08_metaeuk && 
            metaeuk easy-predict ${genfile/${genfile##*.}/softmask.${genfile##*.}} /To/Your/Directory/db/UniRef90/uniref90.metaeukdb $gen.mtek.ur90 /To/Your/Directory/tmp --threads $NP --remove-tmp-files 1 --num-iterations 1 --exhaustive-search 1 --mask-lower-case 1 >> $DIR_out/maker2pred.6.3.log 2>&1 && 
            awk '{OFS="\t"; print $1,$3,$2,$4,$5,$6,$7,$8,$9}' $DIR_out/x08_metaeuk/$gen.mtek.ur90.gff | sed -e "s#Target_ID=#Target=#g;s#TCS_ID=#ID=#g" > $DIR_out/x08_metaeuk/$gen.mtek.ur90.reform.gff && 
            mv $gen.mtek.ur90.fas $gen.mtek.ur90.pep.faa && cp -f $gen.mtek.ur90.pep.faa $DIR_data/busco/pep && 
            mv $gen.mtek.ur90.codon.fas $gen.mtek.ur90.cds.fna && # cp -f $gen.mtek.ur90.cds.fna $DIR_data/busco/cds && 
            echo -e "#######################\nMetaEuk (#1) prediction is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
        else
            echo -e "#######################\nThe file $DIR_out/x08_metaeuk/$gen.mtek.ur90.pep.faa is already there.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
        fi

        cp $DIR_current/maker*.ctl $DIR_out &&
        cd $DIR_out && 

        echo -e "#######################\nChange configuration file for MAKER (run #3) start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
        sed -i -r "s#^genome=#&$genfile#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^maker_gff=#&$DIR_out\/x05_Run2.maker.output\/Run2.all.gff#g" $DIR_out/maker_opts.ctl && 
        # sed -i -e "s#^est=#&$est_evid#g" $DIR_out/maker_opts.ctl && 
        # sed -i -e "s#^protein=#&$prot_evid#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^est_pass=0#est_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^altest_pass=0#altest_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^protein_pass=0#protein_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^rm_pass=0#rm_pass=1#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^model_pass=0#model_pass=1#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^pred_pass=0#pred_pass=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^other_pass=0#other_pass=1#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^est2genome=1#est2genome=0#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^protein2genome=1#protein2genome=0#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^snaphmm=#&$DIR_out\/x06_snap2\/$gen.snap2.hmm#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^gmhmm=#&$DIR_out\/genemarkes2\/$gen.gm2.hmm#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^augustus_species=#&$gen.Run2:self,ancylostoma_ceylanicum:ancylostoma,brugia:brugia,caenorhabditis:caenorhabditis,c_elegans_trsk:c_elegans,trichinella:trichinella#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^pred_gff=#&$DIR_out\/x08_metaeuk\/$gen.mtek.ur90.reform.gff#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^run_evm=0#run_evm=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^pred_stats=0#pred_stats=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^alt_splice=0#alt_splice=1#g" $DIR_out/maker_opts.ctl && 
        sed -i -r "s#^keep_preds=0#keep_preds=1#g" $DIR_out/maker_opts.ctl && 
        # sed -i -r "s#^correct_est_fusion=1#correct_est_fusion=0#g" $DIR_out/maker_opts.ctl && 
        echo -e "#######################\nChange configuration file for MAKER (run #3) processes are done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 

        echo -e "#######################\nMAKER (run #3) start...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 
        /sw/lib/openmpi-4.0.2-gcc-4.8.5/bin/mpiexec -np $NP --oversubscribe maker --ignore_nfs_tmp -base Run3 -RM_off >> $DIR_out/maker2pred.6.3.log 2>&1 && 
        mv $DIR_out/Run3.maker.output $DIR_out/x09_Run3.maker.output && 
        echo -e "#######################\nMAKER (run #3) processes are done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 

        if [[ $(ls -l $DIR_out/x09_Run3.maker.output/Run3_master_datastore_index.log | awk '{ print $5 }') -gt 0 ]]; then

            cd $DIR_out/x09_Run3.maker.output && 

            echo -e "#######################\nCombining MAKER (run #3) results...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log && 

            gff3_merge -d $DIR_out/x09_Run3.maker.output/Run3_master_datastore_index.log >> $DIR_out/maker2pred.6.3.log 2>&1 && 
            fasta_merge -d $DIR_out/x09_Run3.maker.output/Run3_master_datastore_index.log >> $DIR_out/maker2pred.6.3.log 2>&1 && 
            # get rid of fasta sequences in the gff file
            fastastart=$(($(sed -n '/##FASTA/=' Run3.all.gff)-1)) && 
            sed -n "1,$fastastart"p Run3.all.gff > Run3.all.nofasta.gff && 
            maker_map_ids --prefix maker_ --justify 8 --iterate 1 Run3.all.nofasta.gff > Run3.all.idmap && 
            cp Run3.all.gff ShortId.Run3.all.gff &&             
            map_gff_ids Run3.all.idmap ShortId.Run3.all.gff && 
            cp Run3.all.nofasta.gff ShortId.Run3.all.nofasta.gff && 
            map_gff_ids Run3.all.idmap ShortId.Run3.all.nofasta.gff && 
            cp Run3.all.maker.proteins.fasta ShortId.Run3.all.maker.proteins.fasta && 
            map_fasta_ids Run3.all.idmap ShortId.Run3.all.maker.proteins.fasta && 
            cp Run3.all.maker.transcripts.fasta ShortId.Run3.all.maker.transcripts.fasta && 
            map_fasta_ids Run3.all.idmap ShortId.Run3.all.maker.transcripts.fasta && 

            rm -rf $DIR_out/x09_Run3.maker.output/mpi_blastdb && 
            rm -rf $DIR_out/x09_Run3.maker.output/Run3_datastore && 
            rm -f $DIR_out/x09_Run3.maker.output/.NFSLock.* 

            echo -e "#######################\nMAKER (run #3) job is done.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
        else
            echo -e "#######################\nMAKER (run #3) job is finished but with no gain.\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log  
        fi
    else
        echo -e "#######################\nNo Run2.all.gff found...\n$(date)\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
    fi
fi

# source $DIR_cnd/bin/deactivate

if [ -f $DIR_out/x09_Run3.maker.output/Run3.all.maker.proteins.fasta ] && [ -f $DIR_out/x09_Run3.maker.output/Run3.all.maker.transcripts.fasta ]; then 
    cp $DIR_out/x09_Run3.maker.output/Run3.all.maker.proteins.fasta $DIR_out/$gen.final.maker.pep.faa && 
    cp $DIR_out/x09_Run3.maker.output/Run3.all.maker.transcripts.fasta $DIR_out/$gen.final.maker.cds.fna && 
    cp $DIR_out/x09_Run3.maker.output/Run3.all.maker.augustus.proteins.fasta $DIR_out/$gen.final.mk-augu.pep.faa && 
    cp $DIR_out/x09_Run3.maker.output/Run3.all.maker.augustus.transcripts.fasta $DIR_out/$gen.final.mk-augu.cds.fna && 
    
    mkdir -p $DIR_data/busco/pep 
    mkdir -p $DIR_data/busco/cds 
    cp $DIR_out/$gen.final.maker.pep.faa $DIR_data/busco/pep && 
    cp $DIR_out/$gen.final.maker.cds.fna $DIR_data/busco/cds && 
    cp $DIR_out/$gen.final.mk-augu.pep.faa $DIR_data/busco/pep && 
    cp $DIR_out/$gen.final.mk-augu.cds.fna $DIR_data/busco/cds && 

    time1=$(date) && end_sec=$(date --date="$time1" +%s) && 
    dsec=$((end_sec-start_sec)) && dhour=$(echo "scale=4; $dsec / 3600" | bc) && 
    echo -e "#######################\nFinal protein and CDS sequences copied to the BUSCO directory, MAKER annotation accomplished.\n$time1\nOverall elapsed time: $dsec sec (= $dhour hr).\n#######################\n\n" >> $DIR_out/maker2pred.6.3.log 
fi
