#!/bin/bash
# Created by the VLSCI job script generator for SLURM on x86
# Mon Feb 09 2015 22:01:44 GMT+1100 (AEDT)

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Partition for the job:
#SBATCH -p main

# The name of the job:
#SBATCH --job-name="hitime-cpp"

# Maximum number of CPU cores used by the job:
#SBATCH --ntasks=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem-per-cpu=24576

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-24:0:00

# Run the job from the directory where it was launched (default):
# The job command(s):

module load boost-gcc/1.57.0
module load pwiz-gcc/3.0.7069
#./hitime.out data/testing.mzML testing.txt
./hitime.out data/rat.mzML rat.out
