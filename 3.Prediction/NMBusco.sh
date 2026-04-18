#!/bin/bash
#SBATCH --job-name=NMBU.x
#SBATCH --exclusive
#SBATCH --partition=mwvdk
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=350GB

NP=$SLURM_CPUS_PER_TASK # 8
NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE] # 16

# scores_cutoff_eff=0.8
lengths_cutoff_eff=2

cd /To/Your/Directory/nematoda/data/busco

DIR_current=$(pwd)
DIR_cnd=/To/Your/Directory/anaconda3
DIR_db=$DIR_cnd/envs/buscoenv/share/buscodb/lineages/nematoda_odb10 
UniteBUSCOs=/To/Your/Directory/nematoda/data/UniteBUSCOs.py 
chmod -R 777 $UniteBUSCOs 
# export PATH=/To/Your/Directory/anaconda3/bin:/To/Your/Directory/anaconda3/envs/buscoenv/bin:$PATH
# export LIBDIR=/To/Your/Directory/anaconda3/envs/buscoenv/share/buscodb/lineages
mkdir -p $DIR_current/pep
mkdir -p $DIR_current/res
mkdir -p $DIR_current/log
mkdir -p $DIR_current/tmp
mkdir -p $DIR_current/united

echo -e "#######################\nBUSCO extraction start...\n$(date)\n#######################\n\n" >> $DIR_current/busco.log && 
echo -e "JOB NAME: $SLURM_JOB_NAME\nJob ID: $SLURM_JOBID\nAllocate Nodes: $SLURM_JOB_NODELIST\n" >> $DIR_current/busco.log && 
time0=$(date) && start_sec=$(date --date="$time0" +%s) && 
echo -e "#######################\nConfiguring environment...\n$(date)\n#######################\n\n" >> $DIR_current/busco.log 

summarizeBUSCO(){
    local gen=$1 
    shortSummary=$2 
    if [ ! -f $3 ]; then
        echo "Genome,ExecDate,TotalPCGs,Database,DBUpdate,NumGenomes,TotalBUSCOs,nComp,%Comp,nSing,%Sing,nDupl,%Dupl,nFrag,%Frag,nMiss,%Miss" > $3 
    fi
    l2=$(grep -n "The lineage dataset is:" $shortSummary)
    ed=$(date +%Y-%m-%d)
    db=$(echo $l2 | awk '{ print $6 }') 
    dt=$(echo $l2 | awk '{ print $9 }' | sed "s/,//g") 
    ngm=$(echo $l2 | awk '{ print $13 }' | sed "s/,//g") 
    totBusco=$(echo $l2 | awk '{ print $17 }' | sed "s/)//g") 
    l8=$(sed -n $(($(grep -n "***** Results: *****" $shortSummary | sed -e "s/:/ /g" | awk '{ print $1 }') + 2))p $shortSummary | sed -e "s/C://g;s/%\[S:/ /g;s/%,D:/ /g;s/%\],F:/ /g;s/%,M:/ /g;s/%,n:/ /g") 
    pC=$(echo $l8 | awk '{ print $1 }') 
    pS=$(echo $l8 | awk '{ print $2 }') 
    pD=$(echo $l8 | awk '{ print $3 }') 
    pF=$(echo $l8 | awk '{ print $4 }') 
    pM=$(echo $l8 | awk '{ print $5 }') 
    nC=$(grep -n "Complete BUSCOs (C)" $shortSummary | awk '{ print $2 }') 
    nS=$(grep -n "Complete and single-copy BUSCOs (S)" $shortSummary | awk '{ print $2 }') 
    nD=$(grep -n "Complete and duplicated BUSCOs (D)" $shortSummary | awk '{ print $2 }') 
    nF=$(grep -n "Fragmented BUSCOs (F)" $shortSummary | awk '{ print $2 }') 
    nM=$(grep -n "Missing BUSCOs (M)" $shortSummary | awk '{ print $2 }') 
    echo -e "$gen,$ed,$4,$db,$dt,$ngm,$totBusco,$nC,$pC,$nS,$pS,$nD,$pD,$nF,$pF,$nM,$pM" >> $3 
}

echo -e "=======================\nChange OrthoDB configuration...\n$(date)\nscores_cutoff: * $scores_cutoff_eff\nlengths_cutoff: * $lengths_cutoff_eff\n=======================\n\n" >> $DIR_current/busco.log 
# rm $DIR_db/scores_cutoff && 
rm $DIR_db/lengths_cutoff && 

# awk -v x=$scores_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; scr = $2; nscr = scr * x; print gen, nscr }' $DIR_db/scores_cutoff.ori > $DIR_db/scores_cutoff && # used with caution in that it may include more paralogs.
awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; B = $2; delta = $3; len = $4; ndelta = delta * x; print gen, B, ndelta, len }' $DIR_db/lengths_cutoff.ori > $DIR_db/lengths_cutoff # if you use this line, the lengths_cutoff_eff should be greater than 1, recommended as 2; this method was used by BUSCO; but the two are almost equal in the effect.
# awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; B = $2; delta = $3; len = $4; nlen = len * x; print gen, B, delta, nlen }' $DIR_db/lengths_cutoff.ori > $DIR_db/lengths_cutoff # if you use this line, the lengths_cutoff_eff should be in [0,1], recommended as 0.8.

source $DIR_cnd/bin/activate $DIR_cnd/envs/buscoenv && echo "BUSCO environment activated." >> $DIR_current/busco.log || echo "Wrong with BUSCO activation." >> $DIR_current/busco.log 

# Reading input data list
genofile=$(ls $DIR_current/pep/*.pep.faa | sed '/final.mk-augu/d')

if [[ -f $DIR_current/NematodaBuscoDoneList.dat ]]; then
    donefile0=$(cat $DIR_current/NematodaBuscoDoneList.dat) 
    # rm $DIR_current/NematodaBuscoDoneList.dat
else
    donefile0="" 
fi

# MultiTask preparation (set "i<=N" to allow N tasks (a loop refers to a task) at same time and let others waiting in a queue)
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile
for ((i=1;i<=$NT;i++))
do
    echo >&9
done

# BUSCO block
for g in $genofile
do 
    read -u9
    {

    gen=$(echo $g | sed -e "s#.pep.faa##g;s#$DIR_current\/pep\/##g") && 
    totPCGs=$(grep '>' $g | wc -l) && 
    DIR_out=$DIR_current/tmp/$gen && 
    if [[ ! $donefile0 =~ "$gen FINISHED" ]]; then
        echo -e "$(date) $gen START" >> $DIR_current/NematodaBuscoDoneList.dat &&         
        echo -e "#######################\nStart BUSCO for $gen\nTotally $totPCGs proteins.\n$(date)\n#######################\n\n" >> $DIR_current/busco.log 
        # source $DIR_cnd/bin/activate $DIR_cnd/envs/buscoenv && 
        # echo "BUSCO environment activated." >> $DIR_current/busco.log && 
        cd $DIR_current/tmp && 
        busco -m protein --cpu $NP --offline -i $g -o $gen -l $DIR_db >> $DIR_current/busco.log && 
        summarizeBUSCO $gen $DIR_out/short_summary.specific.nematoda_odb10.$gen.txt $DIR_current/allBUSCOsummary.csv $totPCGs && 
        cp -r $DIR_out/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences $DIR_current/res/$gen && 
        cp -r $DIR_out/logs $DIR_current/log/$gen && 
        echo -e "$(date) $gen FINISHED" >> $DIR_current/NematodaBuscoDoneList.dat && 
        rm -rf $DIR_out 
    fi 

    echo >&9 
    } & 
done

wait 

donefile1=$(cat $DIR_current/NematodaBuscoDoneList.dat) 
taskstarted=$(grep -n 'START' $DIR_current/NematodaBuscoDoneList.dat | awk '{ print $7 }') 
sjob=$(grep -n 'START' $DIR_current/NematodaBuscoDoneList.dat | wc -l)

fjob=0
for j in $taskstarted
do 
    if [[ $donefile1 =~ "$j FINISHED" ]]; then
        let fjob++ 
    else
        echo -e "$j" >> $DIR_current/UnfinishedJob.dat
    fi
done

echo -e "=======================\nRecover OrthoDB default configuration...\n$(date)\n=======================\n\n" >> $DIR_current/busco.log 
cp -f $DIR_db/scores_cutoff.ori $DIR_db/scores_cutoff && 
cp -f $DIR_db/lengths_cutoff.ori $DIR_db/lengths_cutoff && 

ujob=$((sjob-fjob)) && 

time1=$(date) && end_sec=$(date --date="$time1" +%s) && 
dsec=$((end_sec-start_sec)) && dhour=$(echo "scale=4; $dsec / 3600" | bc) 

if [ $fjob -eq $sjob ]; then
    echo -e "#######################\nAll BUSCO extraction jobs are accomplished.\n$time1\nOverall elapsed time: $dsec sec (= $dhour hr).\n#######################\n\n" >> $DIR_current/busco.log 
else
    echo -e "#######################\nAll BUSCO jobs are expired, $fjob are finished, $ujob are unfinished. Please check the file $DIR_current/UnfinishedJob.dat for details.\n$time1\nOverall elapsed time: $dsec sec (= $dhour hr).\n#######################\n\n" >> $DIR_current/busco.log 
fi

echo -e "#######################\nStarting to unite the BUSCOs of same species.\n$(date)\n#######################\n\n" >> $DIR_current/busco.log 
cd $DIR_current && 
$UniteBUSCOs --indir $DIR_current/res --outdir $DIR_current/united --sortmissing >> $DIR_current/busco.log 2>&1 && 
echo -e "#######################\nBUSCOs united.\n$(date)\n#######################\n\n" >> $DIR_current/busco.log 

cd $DIR_current
zip -r buscores.zip $DIR_current/res && 
zip -r buscolog.zip $DIR_current/log 
