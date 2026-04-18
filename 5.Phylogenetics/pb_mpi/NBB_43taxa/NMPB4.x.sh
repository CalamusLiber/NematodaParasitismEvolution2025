#!/bin/bash
#SBATCH --job-name=NMPB4
#SBATCH --nodes=3
#SBATCH --ntasks=72
#SBATCH --mem=200GB
#SBATCH --exclusive

DIR_cnd=/To/Your/Directory/anaconda3
#NP=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE*$SLURM_CPUS_PER_TASK]
source $DIR_cnd/bin/activate $DIR_cnd/envs/phylobayes

cd /public1/home/scb4616/AllBackbone4.nom.rmsp

mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi -d ./NBBData43taxa174426AA.bmge.phy -gtr -cat -dgam 4 -x 1 2000 AllBackbone4.nom.rmsp.pb.chain$1  

# # for the second, third, or more 2000 generations, please use the following code:
# mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi AllBackbone4.nom.rmsp.pb.chain$1 

# # for the gtr model without the CAT model, please use the following code:
# mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi -d ./NBBData43taxa174426AA.bmge.phy -gtr -ncat 1 -dgam 4 -x 1 2000 AllBackbone4.nom.rmsp.pb.gtr.chain$1  
# mpirun -np $SLURM_NTASKS --mca btl tcp,self --mca btl_tcp_if_include ib0 pb_mpi AllBackbone4.nom.rmsp.pb.gtr.chain$1 
