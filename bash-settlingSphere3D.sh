#!/bin/bash
#SBATCH --partition=dev_single
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ujvow@student.kit.edu
#SBATCH --job-name=olbench

module purge
module load compiler/gnu/10.2
module load mpi/openmpi/4.1
make clean
make

export WORKSPACE=/pfs/work7/workspace/scratch/ujvow-particleSimulation-0

if [ -d tmp ]; then
  if [ -L tmp ]; then
    unlink tmp
  else
    rm -rf tmp
  fi
fi

mkdir $WORKSPACE/job_bwUC2_${SLURM_JOB_ID}
ln -s $WORKSPACE/job_bwUC2_${SLURM_JOB_ID} tmp

mpirun -n 40 --bind-to core --map-by core ./settlingSphere3D | tee -a settlingSphere.txt

# prints "size, steps, cuboids_per_process, process count, thread count, MLUPs per process, total MLUPs"
#for size in 100. 150. 200. 250.;
#do
#   for n in 1 5 10 20 30 40;
#      do
#              mpirun -n $n --bind-to core --map-by core ./cavity3d $size
#      done
#done

