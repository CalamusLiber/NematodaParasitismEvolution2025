# Shell function for hard-masking FASTA sequences using RepeatMasker-like GFF coordinates.
#
# Inputs:
#   $1: input FASTA file
#   $2: GFF-like annotation file (contig_id ... start end ...)
#
# Output:
#   A new FASTA file where masked intervals are replaced by 'N'.
#   Output name is derived from input, e.g. genome.fna -> genome.hardmask.fna
GffHardMasker(){
    local infna=$1 
    local gff=$2 
    local outfna=${infna/${infna##*.}/hardmask.${infna##*.}} 
    
    # Record FASTA header line indices and contig IDs.
    local LineNum=($(grep -n "^>" $infna | awk -F ":>" '{print $1}')) 
    local ContigTagQ=($(grep "^>" $infna | awk -F "[> ]" '{print $2}')) 
    
    # Process each contig separately, convert sequence to a base-level array,
    # then replace masked positions with 'N'.
    for g in ${!ContigTagQ[*]} 
    do 
        # Extract the sequence of contig g and split it into a character array.
        local sline=$[${LineNum[$g]}+1] 
        local eline=$[${LineNum[$[$g+1]]}-1] 
        local seqx=($(sed -n "${sline},${eline}p" $infna | sed ':t;N;s/\n//;b t' | sed -e "s/./& /g")) 
        
        # Collect masking intervals (1-based) for this contig from the GFF file.
        local Qstart=($(grep "^${ContigTagQ[$g]}" $gff | awk '{ print $4 }')) 
        local Qend=($(grep "^${ContigTagQ[$g]}" $gff | awk '{ print $5 }')) 

        local xpos=()
        for i in ${!Qstart[*]}
        do 
            # Convert each interval to 0-based coordinates for array indexing.
            xpos=(${xpos[*]} $(seq $[${Qstart[$i]}-1] $[${Qend[$i]}-1])) 
        done
        
        # Replace all masked positions with 'N'.
        for j in ${xpos[*]} 
        do 
            seqx[$j]=N 
        done

        seqx=$(echo "${seqx[*]}" | sed "s/ //g") 

        # Write masked sequence to output FASTA (80 bp per line).
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
