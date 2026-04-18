#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --exclusive

cd /To/Your/Directory/nematoda/data/rnaseq
DIR_cur=$(pwd)
DIR_cnd=/To/Your/Directory/anaconda3

source $DIR_cnd/bin/activate $DIR_cnd/envs/rnaseqenv && 

mkdir -p $DIR_cur/Enoplus_brevis_NEB.A64176.soap_trans && 

cd $DIR_cur/Enoplus_brevis_NEB.A64176.soap_trans && 

SOAPdenovo-Trans-127mer all -K $1 -p 24 -s $DIR_cur/soapdenovo_trans.cfg -o res.k$1 >> NMSoapTrans.1.k$1.log 2>&1 

