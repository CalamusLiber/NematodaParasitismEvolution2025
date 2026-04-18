#!/bin/bash
#SBATCH --job-name=NMPBmod
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=360GB

DIR=$(pwd)
NP=64
# DIR_cnd=/public1/home/scb4616/apps/anaconda3
# DIR_cnd=/To/Your/Directory/anaconda3
DIR_cnd=/usr/local/anaconda3 
PY_loocv=$DIR_cnd/pkgs/phylobayes-mpi-1.9-h5c6ebe3_0/scripts/read_loocv_waic.py 
PY=$DIR_cnd/bin/python3.9 
source $DIR_cnd/bin/activate phylobayes 
conda activate phylobayes

cd $DIR/AllBackbone1.nom.rmsp

nohup mpirun -np 48 pb_mpi -d FBBData58taxa149581AA.bmge.phy -T bpcomp1catgtr.con.tre -cat -gtr -dgam 4 -x 1 1200 CatGtrG4.chain1 & 
nohup mpirun -np 48 pb_mpi -d FBBData58taxa149581AA.bmge.phy -T bpcomp1catgtr.con.tre -cat -gtr -dgam 4 -x 1 1200 CatGtrG4.chain2 & 
nohup mpirun -np 24 pb_mpi -d FBBData58taxa149581AA.bmge.phy -T bpcomp1gtr.con.tre -gtr -ncat 1 -dgam 4 -x 1 1200 GtrG4.chain1 & 
nohup mpirun -np 24 pb_mpi -d FBBData58taxa149581AA.bmge.phy -T bpcomp1gtr.con.tre -gtr -ncat 1 -dgam 4 -x 1 1200 GtrG4.chain2 & 

cd $DIR/AllBackbone4.nom.rmsp

nohup mpirun -np 48 pb_mpi -d NBBData43taxa174426AA.bmge.phy -T bpcomp4catgtr.con.tre -cat -gtr -dgam 4 -x 1 1200 CatGtrG4.chain1 & 
nohup mpirun -np 48 pb_mpi -d NBBData43taxa174426AA.bmge.phy -T bpcomp4catgtr.con.tre -cat -gtr -dgam 4 -x 1 1200 CatGtrG4.chain2 & 
nohup mpirun -np 24 pb_mpi -d NBBData43taxa174426AA.bmge.phy -T bpcomp4gtr.con.tre -gtr -ncat 1 -dgam 4 -x 1 1200 GtrG4.chain1 & 
nohup mpirun -np 24 pb_mpi -d NBBData43taxa174426AA.bmge.phy -T bpcomp4gtr.con.tre -gtr -ncat 1 -dgam 4 -x 1 1200 GtrG4.chain2 & 

# Leave-one-out cross-validation and the wAIC
cd $DIR/AllBackbone1.nom.rmsp

mpirun -np 30 readpb_mpi -x 300 9 -sitelogl CatGtrG4.chain1 & 
mpirun -np 30 readpb_mpi -x 300 9 -sitelogl CatGtrG4.chain2 & 
mpirun -np 30 readpb_mpi -x 300 9 -sitelogl GtrG4.chain1 & 
mpirun -np 30 readpb_mpi -x 300 9 -sitelogl GtrG4.chain2 & 

$PY $PY_loocv CatGtrG4.chain?.sitelogl > AllBackbone1.nom.rmsp.pb.CatGtrG4.loocv 2>&1 
$PY $PY_loocv GtrG4.chain?.sitelogl > AllBackbone1.nom.rmsp.pb.gtr.GtrG4.loocv 2>&1 

cd $DIR/AllBackbone4.nom.rmsp

mpirun -np 30 readpb_mpi -x 300 9 -sitelogl CatGtrG4.chain1 & 
mpirun -np 30 readpb_mpi -x 300 9 -sitelogl CatGtrG4.chain2 & 
mpirun -np 30 readpb_mpi -x 300 9 -sitelogl GtrG4.chain1 & 
mpirun -np 30 readpb_mpi -x 300 9 -sitelogl GtrG4.chain2 & 

$PY $PY_loocv CatGtrG4.chain?.sitelogl > AllBackbone4.nom.rmsp.pb.CatGtrG4.loocv 2>&1 
$PY $PY_loocv GtrG4.chain?.sitelogl > AllBackbone4.nom.rmsp.pb.gtr.GtrG4.loocv 2>&1 

# Posterior predictive model checking for adequacy (MCMC chains should be at stationarity)
cd $DIR/AllBackbone1.nom.rmsp

nohup mpirun -np 60 readpb_mpi -x 300 9 -allppred CatGtrG4.chain1 & 
nohup mpirun -np 60 readpb_mpi -x 300 9 -allppred GtrG4.chain1 & 

cd $DIR/AllBackbone4.nom.rmsp

nohup mpirun -np 60 readpb_mpi -x 300 9 -allppred CatGtrG4.chain1 & 
nohup mpirun -np 60 readpb_mpi -x 300 9 -allppred GtrG4.chain1 & 

