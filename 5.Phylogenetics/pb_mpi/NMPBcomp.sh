#!/bin/sh
DIR_cnd=/public1/home/scb4616/apps/anaconda3
#DIR_cnd=/To/Your/Directory/anaconda3
source $DIR_cnd/bin/activate phylobayes

filehead=$(pwd | awk -v FS='/' '{print $5}')
#filehead=$(pwd | awk -v FS='/' '{print $7}')

for run1 in `seq $1 $2`
do 
    if [ ! -f ${filehead}.pb.chain${run1}.trace ]; then 
        continue
    fi 
	
	bpcomp -x $3 ${filehead}.pb.chain${run1} -o bpcomp${run1} > bpcomp${run1}.log 2>&1 
    tracecomp -x $3 ${filehead}.pb.chain${run1} -o tracecomp${run1} > tracecomp${run1}.log 2>&1 
    
	for run2 in `seq $1 $2` 
    do 
        if [ ! -f ${filehead}.pb.chain${run2}.trace ]; then 
            continue
        fi 

        if [ ${run1} -lt ${run2} ]; then 
            bpcomp -x $3 ${filehead}.pb.chain${run1} ${filehead}.pb.chain${run2} -o bpcomp${run1}+${run2} > bpcomp${run1}+${run2}.log 2>&1 
            tracecomp -x $3 ${filehead}.pb.chain${run1} ${filehead}.pb.chain${run2} -o tracecomp${run1}+${run2} > tracecomp${run1}+${run2}.log 2>&1 
        
			for run3 in `seq $1 $2` 
			do 
				if [ ! -f ${filehead}.pb.chain${run3}.trace ]; then 
					continue
				fi 
			
				if [ ${run2} -lt ${run3} ]; then 
					bpcomp -x $3 ${filehead}.pb.chain${run1} ${filehead}.pb.chain${run2} ${filehead}.pb.chain${run3} -o bpcomp${run1}+${run2}+${run3} > bpcomp${run1}+${run2}+${run3}.log 2>&1 
					tracecomp -x $3 ${filehead}.pb.chain${run1} ${filehead}.pb.chain${run2} ${filehead}.pb.chain${run3} -o tracecomp${run1}+${run2}+${run3} > tracecomp${run1}+${run2}+${run3}.log 2>&1 
				
					for run4 in `seq $1 $2` 
					do 
						if [ ! -f ${filehead}.pb.chain${run4}.trace ]; then 
							continue
						fi 
			
						if [ ${run3} -lt ${run4} ]; then 
							bpcomp -x $3 ${filehead}.pb.chain${run1} ${filehead}.pb.chain${run2} ${filehead}.pb.chain${run3} ${filehead}.pb.chain${run4} -o bpcomp${run1}+${run2}+${run3}+${run4} > bpcomp${run1}+${run2}+${run3}+${run4}.log 2>&1 
							tracecomp -x $3 ${filehead}.pb.chain${run1} ${filehead}.pb.chain${run2} ${filehead}.pb.chain${run3} ${filehead}.pb.chain${run4} -o tracecomp${run1}+${run2}+${run3}+${run4} > tracecomp${run1}+${run2}+${run3}+${run4}.log 2>&1 
						fi 
					done
				fi 
			done
		fi 
    done 
done 

grep 'maxdiff' bpcomp*.log > Maxdiff_report.txt 
