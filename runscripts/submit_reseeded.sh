#!/bin/bash
#SBATCH --job-name b2.4_1.234B
#SBATCH --qos standard
#SBATCH --time 2-0:00:00
#SBATCH --account dp208
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --dependency=singleton

module purge
module load cuda/11.4.1  openmpi/4.1.1-cuda11.4.1  ucx/1.12.0-cuda11.4.1

export OMP_NUM_THREADS=4
export OMPI_MCA_btl=^uct,openib
export UCX_TLS=gdr_copy,rc,rc_x,sm,cuda_copy,cuda_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_IB_GPU_DIRECT_RDMA=yes
export UCX_MEMTYPE_CACHE=n

VOL=32.32.32.64
MPI=2.2.1.2
TRAJECTORIES=500
BETA=2.4
MINUSMASS=1.234
NSTEPS=50
TLEN=2.0
SAVEFREQ=1
STARTTRAJ=$(ls -rt cnfg/ckpoint_lat.*[^k] | tail -1 | sed 's/[^0-9]*//')

mpirun -np ${SLURM_NTASKS} -x LD_LIBRARY_PATH  \
       --bind-to none ./wrappersimple.sh ~/src/Grid/build/tests/hmc/Test_rhmc_Wilson1p1 \
       --StartingType CheckpointStart \
       --StartingTrajectory ${STARTTRAJ} \
       --grid ${VOL} \
       --mpi ${MPI} \
       --accelerator-threads 8 \
       --Trajectories ${TRAJECTORIES} \
       --Thermalizations 0 \
       --beta ${BETA} \
       --fermionmass -${MINUSMASS} \
       --nsteps ${NSTEPS} \
       --tlen ${TLEN} \
       --savefreq ${SAVEFREQ} \
       --serialseed "11 12 13 14 15" \
       --parallelseed "16 17 18 19 20" > hmc_${SLURM_JOB_ID}.out

