#!/bin/bash
#----------- Warning: Before executing this script, ---------------------------#
#----------- make sure that there are NO SumMCMCT jobs (SumMCMCtree.py) -----------#
#----------- running in the background or staying in the squeue! --------------#

echo -e "=============\nWarning: Before executing this script, \nmake sure that there are NO SumMCMCT jobs (SumMCMCtree.py) \nrunning in the background or staying in the squeue!  \n=============\n"

DIR_rt=$(pwd)

# Extract topology name from current directory path.
# Example: /path/to/Nema.XYZ.BD -> Topo=XYZ
Topo=$(echo $(pwd) | awk -v FS='/' '{print $NF}' | awk -v FS='.BD' '{print $1}')

# Build the list of calibration names from the .tre files in the current directory.
# Expected file names: Nema.<Topo>.<Calibration>.tre
TreeCalibration=$(ls *.tre | sed "s#\/##g;s#Nema.${Topo}.##g;s#.tre##g")

# Define the clock models and partition schemes used in the analysis.
ClockModel='IR AR'      # IR = independent rates, AR = autocorrelated rates
Partition='1P 5P'       # 1P = one partition, 5P = five-partition strategy
rgene_gamma='15 20 25'  # Gamma rate parameter choices for rgene distribution

# Iterate through each analysis combination and cleanup residual directories.
for TreCali in $TreeCalibration
 do
    for Part in $Partition
    do
        for Clock in $ClockModel
        do
            for Rgene in $rgene_gamma
            do
                DIR_cur=$DIR_rt/${TreCali}.${Part}.${Clock}.rg${Rgene}
                cd $DIR_cur

                # When a residual output directory exists alongside the archived tarball,
                # compress the directory to ZIP, remove the old tarball, and delete the directory.
                if [ -d DIR_DIV -a -f DIR_DIV.tar.gz ]
                then
                    echo "Clean residual DIR_DIV folder in $DIR_cur"
                    zip -r DIR_DIV.zip DIR_DIV
                    rm DIR_DIV.tar.gz
                    rm -r DIR_DIV
                fi

                if [ -d DIR_PRD -a -f DIR_PRD.tar.gz ]
                then
                    echo "Clean residual DIR_PRD folder in $DIR_cur"
                    zip -r DIR_PRD.zip DIR_PRD
                    rm DIR_PRD.tar.gz
                    rm -r DIR_PRD
                fi
            done
        done
    done
 done
