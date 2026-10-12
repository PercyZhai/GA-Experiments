#!/bin/bash
#---------------------------------------------------------------------------------
# simulation/submit.sh
# Generic SLURM job script: runs one cell (n, p, beta, setup) of the simulation
# grid (paper Sec. VII and Appendix I) with simulation/lowdim_sim.R or
# simulation/highdim_sim.R.
#
# Usage (from the repository root):
#   mkdir -p logs
#   sbatch simulation/submit.sh <lowdim|highdim> <n> <p> <beta> [toeplitz|banded]
#   e.g. sbatch simulation/submit.sh lowdim 1000 10 0.55 toeplitz
# simulation/run_grid.sh submits the whole grid.
#
# Outputs: results/lowdim/ or results/highdim/ (written by the R script),
#          logs/<regime>_<n>_<p>_<beta>_<setup>.Rout (R output),
#          logs/%x_%j.out, logs/%x_%j.err (SLURM stdout/stderr; %x = job name, %j = job ID).
# logs/ must exist before submission, since SLURM opens the --output/--error
# files before this script starts.
#
# Adjust the #SBATCH lines below to your cluster. The high-dimensional cells run in
# parallel on SLURM_CPUS_PER_TASK cores and need more CPUs, memory and time (see
# simulation/highdim_sim.R); options given to sbatch override the lines below, e.g.
#   sbatch --cpus-per-task=28 --mem=200G --time=3-00:00:00 simulation/submit.sh highdim 800 8000 0.55
# (simulation/run_grid.sh passes such options for the high-dimensional cells).
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Account information (uncomment and set if your cluster requires an account)

##SBATCH --account=<your_account>

#---------------------------------------------------------------------------------
# Resources requested

##SBATCH --partition=<your_partition>   # uncomment and set if your cluster requires a partition
#SBATCH --cpus-per-task=1          # number of CPUs requested
#SBATCH --mem=32G                  # requested memory (per job)
#SBATCH --time=1-00:00:00          # wall clock limit (d-hh:mm:ss)

#---------------------------------------------------------------------------------
# Job specific name (helps organize and track progress of jobs)

#SBATCH --job-name=ga_sim          # run_grid.sh sets one name per cell

#---------------------------------------------------------------------------------
# Output Management
# Direct slurm system output to the 'logs' folder
# %x = Job Name, %j = Job ID. This prevents overwriting.
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#---------------------------------------------------------------------------------
# Print some useful variables

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

#---------------------------------------------------------------------------------
# Load necessary modules for the job (example only; module names differ between clusters)

# module load R

#---------------------------------------------------------------------------------
# Strict Argument Checking
# Check if 4 arguments are provided. If not, exit with error.
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
  echo "ERROR: Missing arguments."
  echo "Usage: sbatch simulation/submit.sh <lowdim|highdim> <n> <p> <beta> [toeplitz|banded]"
  exit 1
fi

# Assign arguments to variables
REGIME=$1
N_VAL=$2
P_VAL=$3
BETA_VAL=$4
SETUP=${5:-toeplitz}

if [ "$REGIME" != "lowdim" ] && [ "$REGIME" != "highdim" ]; then
  echo "ERROR: regime must be lowdim or highdim, got '$REGIME'."
  exit 1
fi

if [ "$SETUP" != "toeplitz" ] && [ "$SETUP" != "banded" ]; then
  echo "ERROR: setup must be toeplitz or banded, got '$SETUP'."
  exit 1
fi

if [ ! -f functions/bootstrap_sample.R ]; then
  echo "ERROR: Submit this script from the root of the GA-Experiments repository."
  exit 1
fi

echo "Running job with: regime=$REGIME, n=$N_VAL, p=$P_VAL, beta=$BETA_VAL, setup=$SETUP"

#---------------------------------------------------------------------------------
# Execution (in the directory sbatch was called from, i.e. the repository root)

# 1. Create the logs directory if it doesn't exist
mkdir -p logs

# 2. Define log file name dynamically
LOG_FILE="logs/${REGIME}_${N_VAL}_${P_VAL}_${BETA_VAL}_${SETUP}.Rout"

# 3. Run R
Rscript "simulation/${REGIME}_sim.R" "$N_VAL" "$P_VAL" "$BETA_VAL" "$SETUP" > "$LOG_FILE" 2>&1
