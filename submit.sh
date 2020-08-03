#!/bin/bash --login
#SBATCH -D /home/fh04/fh04/bvanbun/model/    # the working directory
#SBATCH -o /home/fh04/fh04/bvanbun/model/output/output.%A.out    # file for stdout
#SBATCH -e /home/fh04/fh04/bvanbun/model/output/output.%A.err    # file for stderr
#SBATCH -J SEEEIRDS_CH    # job name
#SBATCH -p workq          # name of the queue - do not change
#SBATCH --nodes=16         # number of nodes you want to use
#SBATCH --tasks-per-node=64    # number of MPI processes per node (64 is standard)
#SBATCH --cpus-per-task=1      # do not change

# Load any required modules
module load PrgEnv/gnu-9.2/openmpi-4.0.2

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

# -n is the total number of MPI processes - here: 5 x 64
mpirun -n 1000 --mca routed radix bin/al3c config/SEEEIRDS_CH.xml
