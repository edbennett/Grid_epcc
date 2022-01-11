#!/bin/bash
#SBATCH --job-name b2.4_48x24
#SBATCH --qos standard
#SBATCH --time 2-0:00:00
#SBATCH --account dp208
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module purge
module load cuda/11.4.1  openmpi/4.1.1-cuda11.4.1  ucx/1.12.0-cuda11.4.1

export OMP_NUM_THREADS=4
export OMPI_MCA_btl=^uct,openib
export UCX_TLS=gdr_copy,rc,rc_x,sm,cuda_copy,cuda_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_IB_GPU_DIRECT_RDMA=yes
export UCX_MEMTYPE_CACHE=n


cd m1.21_48x24/run1
STARTTRAJ=$(($(ls -rt cnfg/ckpoint_lat.*[^k] | tail -1 | sed 's/[^0-9]*//')))
NSTEPS=45
mpirun -np 1 -x LRANK=0 -x LD_LIBRARY_PATH  --bind-to none ./wrapper1gpu.sh ~/src/Grid/build/tests/hmc/Test_rhmc_Wilson1p1 --StartingType CheckpointStart --StartingTrajectory ${STARTTRAJ} --grid 24.24.24.48 --mpi 1.1.1.1  --accelerator-threads 8 --Trajectories 500 --Thermalizations 0 --beta 2.4 --fermionmass -1.21 --nsteps ${NSTEPS} --tlen 2.0 --savefreq 1 > hmc_${SLURM_JOB_ID}.out &
cd -

cd m1.21_48x24/run2
STARTTRAJ=$(($(ls -rt cnfg/ckpoint_lat.*[^k] | tail -1 | sed 's/[^0-9]*//')))
NSTEPS=45
mpirun -np 1 -x LRANK=1 -x LD_LIBRARY_PATH  --bind-to none ./wrapper1gpu.sh ~/src/Grid/build/tests/hmc/Test_rhmc_Wilson1p1 --StartingType CheckpointStart --StartingTrajectory ${STARTTRAJ} --grid 24.24.24.48 --mpi 1.1.1.1  --accelerator-threads 8 --Trajectories 500 --Thermalizations 0 --beta 2.4 --fermionmass -1.21 --nsteps ${NSTEPS} --tlen 2.0 --savefreq 1 --serialseed "11 12 13 14 15" --parallelseed "16 17 18 19 20" > hmc_${SLURM_JOB_ID}.out &
#mpirun -np 1 -x LRANK=1 -x LD_LIBRARY_PATH  --bind-to none ./wrapper1gpu.sh ~/src/Grid/build/tests/hmc/Test_rhmc_Wilson1p1 --StartingType CheckpointStartWithReseed --StartingTrajectory ${STARTTRAJ} --grid 24.24.24.48 --mpi 1.1.1.1  --accelerator-threads 8 --Trajectories 500 --Thermalizations 0 --beta 2.4 --fermionmass -1.21 --nsteps ${NSTEPS} --tlen 2.0 --savefreq 1 --serialseed "11 12 13 14 15" --parallelseed "16 17 18 19 20" > hmc_${SLURM_JOB_ID}.out &
cd -

cd m1.222_48x24/run1
STARTTRAJ=$(($(ls -rt cnfg/ckpoint_lat.*[^k] | tail -1 | sed 's/[^0-9]*//')))
NSTEPS=50
mpirun -np 1 -x LRANK=2 -x LD_LIBRARY_PATH  --bind-to none ./wrapper1gpu.sh ~/src/Grid/build/tests/hmc/Test_rhmc_Wilson1p1 --StartingType CheckpointStart --StartingTrajectory ${STARTTRAJ} --grid 24.24.24.48 --mpi 1.1.1.1  --accelerator-threads 8 --Trajectories 200 --Thermalizations 0 --beta 2.4 --fermionmass -1.222 --nsteps ${NSTEPS} --tlen 2.0 --savefreq 1 > hmc_${SLURM_JOB_ID}.out &
cd -

cd m1.222_48x24/run2
STARTTRAJ=$(($(ls -rt cnfg/ckpoint_lat.*[^k] | tail -1 | sed 's/[^0-9]*//')))
NSTEPS=50
mpirun -np 1 -x LRANK=3 -x LD_LIBRARY_PATH  --bind-to none ./wrapper1gpu.sh ~/src/Grid/build/tests/hmc/Test_rhmc_Wilson1p1 --StartingType CheckpointStart --StartingTrajectory ${STARTTRAJ} --grid 24.24.24.48 --mpi 1.1.1.1  --accelerator-threads 8 --Trajectories 200 --Thermalizations 0 --beta 2.4 --fermionmass -1.222 --nsteps ${NSTEPS} --tlen 2.0 --savefreq 1 --serialseed "11 12 13 14 15" --parallelseed "16 17 18 19 20" > hmc_${SLURM_JOB_ID}.out &
#mpirun -np 1 -x LRANK=3 -x LD_LIBRARY_PATH  --bind-to none ./wrapper1gpu.sh ~/src/Grid/build/tests/hmc/Test_rhmc_Wilson1p1 --StartingType CheckpointStartWithReseed --StartingTrajectory ${STARTTRAJ} --grid 24.24.24.48 --mpi 1.1.1.1  --accelerator-threads 8 --Trajectories 500 --Thermalizations 0 --beta 2.4 --fermionmass -1.222 --nsteps ${NSTEPS} --tlen 2.0 --savefreq 1 --serialseed "11 12 13 14 15" --parallelseed "16 17 18 19 20" > hmc_${SLURM_JOB_ID}.out &
cd -

wait
