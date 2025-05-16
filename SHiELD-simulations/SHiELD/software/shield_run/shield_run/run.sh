#!/bin/bash
#SBATCH --qos=urgent
#SBATCH --account=gfdl_w
#SBATCH --clusters=c5
#SBATCH --job-name=SHiELD
#SBATCH --output=slurm-logs/%x.%j.out

set -u
set -e
set -v

SCRIPT=$1
EXECUTABLE=$2
SEGMENTS=$3
DESTINATION=$4

# Print the location of the script used
echo "SCRIPT is"
echo ${SCRIPT}

# Activate conda environment
source /sw/gaea-c5/python/3.9/anaconda-base/etc/profile.d/conda.sh
echo "Activating conda environment ${CONDA_DEFAULT_ENV}:"
conda activate ${CONDA_DEFAULT_ENV}

# List environment for simulation
echo "Active conda environment:"
conda list

# Set environment variables following SHiELD runscripts
export MPICH_ENV_DISPLAY
export MPICH_MPIIO_CB_ALIGN=2
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=536870912
export NC_BLKSZ=1M

# necessary for OpenMP when using Intel
export KMP_STACKSIZE=256m
export SLURM_CPU_BIND=verbose

# Required for running on C5 as of 2024-10-16; see
# https://github.com/NOAA-GFDL/SHiELD_build/pull/48
export FI_VERBS_PREFER_XRC=0

# Run the simulation.  Note the -u is required for instantaneous output to
# stdout
python -u -m shield_run.append ${EXECUTABLE} ${DESTINATION}

# Resubmit to run another segment if needed
segment_number=$(python -m shield_run.get_segment_number ${DESTINATION})
if [ ${segment_number} -lt ${SEGMENTS} ]
then
  echo "Resubmitting to run another segment... "
  sbatch \
      --ntasks=${SLURM_NTASKS} \
      --export=SLURM_NTASKS,SBATCH_TIMELIMIT,CONDA_DEFAULT_ENV \
      ${SCRIPT} \
      ${SCRIPT} \
      ${EXECUTABLE} \
      ${SEGMENTS} \
      ${DESTINATION}
  sleep 60
fi
