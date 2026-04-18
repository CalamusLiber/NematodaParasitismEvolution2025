#!/bin/bash
#SBATCH --job-name=NMPB1
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --mem=500GB
#SBATCH --exclusive

DIR_cnd=/public1/home/scb4616/apps/anaconda3
source $DIR_cnd/bin/activate phylobayes

cd /public1/home/scb4616/AllBackbone1.nom.rmsp

mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi -d ./FBBData58taxa149581AA.bmge.phy -gtr -cat -dgam 4 -x 1 2000 AllBackbone1.nom.rmsp.pb.chain$1  

# # for the second, third, or more 2000 generations, please use the following code:
# mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi AllBackbone1.nom.rmsp.pb.chain$1 

# # for the gtr model without the CAT model, please use the following code:
# mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi -d ./FBBData58taxa149581AA.bmge.phy -gtr -ncat 1 -dgam 4 -x 1 2000 AllBackbone1.nom.rmsp.pb.gtr.chain$1  
# mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi AllBackbone1.nom.rmsp.pb.gtr.chain$1 
