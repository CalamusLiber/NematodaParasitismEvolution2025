#!/bin/sh

# soap2trinity.sh
# Convert SOAP assembler contig output into a Trinity-style FASTA file
# Usage: soap2trinity.sh <soap_contig_file> [output_fasta]
# If the output path is omitted, the script writes to the same folder as the input

soapin=$1

# Determine output file path
if [ "$2" = "" ]; then 
    if [[ $soapin =~ '/' ]]; then 
        trinityout=${soapin/.soap/.soap_trinity.Trinity}
    else
        trinityout="./${soapin/.soap/.soap_trinity.Trinity}"
    fi
else
    if [[ $2 =~ '/' ]]; then 
        trinityout=$2
    else
        trinityout="./$2"
    fi
fi

# Find line numbers of SOAP contig headers (lines beginning with '>')
LineNum=($(grep -n "^>" $soapin | awk -F ":>" '{print $1}')) 
# Convert SOAP contig IDs to Trinity-compatible contig names
ContigTagQ=($(sed -n -e 's/scaffold[0-9]* Locus_/SOAP_DN/g;s/_\([0-9]*\) /_c0_g1_i\1 /g;s/  / /;s/>//p' $soapin | awk '{print $1}')) 
# Total number of lines in the SOAP input file
LastLine=$(cat $soapin | wc -l)

# Convert each SOAP contig into a single-line Trinity-style FASTA entry
for ((g=1;g<${#ContigTagQ[*]};g++)) 
do 
    # Extract the sequence block between header lines for the current contig
    sline=$[${LineNum[$[$g-1]]}+1] 
    eline=$[${LineNum[$[$g]]}-1] 
    seq0=$(sed -n "${sline},${eline}p" $soapin) 
    # Remove whitespace so each contig sequence is one continuous line
    seq1=$(echo $seq0 | sed -e "s/ //g") 

    # Write the contig header and sequence to the output FASTA
    echo ">${ContigTagQ[$g]} len=${#seq1} path=[0:0-$[${#seq1}-1]]" >> $trinityout 
    echo ${seq1} >> $trinityout 
done

# Process the final contig separately using the last header position
seq0=$(sed -n "$[${LineNum[-1]}+1],${LastLine}p" $soapin) 
seq1=$(echo $seq0 | sed -e "s/ //g") 
echo ">${ContigTagQ[-1]} len=${#seq1} path=[0:0-$[${#seq1}-1]]" >> $trinityout 
echo ${seq1} >> $trinityout 
