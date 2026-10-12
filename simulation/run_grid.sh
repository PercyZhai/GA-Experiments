#!/bin/bash
#---------------------------------------------------------------------------------
# simulation/run_grid.sh
# Runs the full simulation grid of the paper (Sec. VII: toeplitz; Appendix I: banded).
#   low-dim:  n in {200, 500, 1000, 2000, 4000, 10000} x p in {2, 10, 30, 100}
#   high-dim: n in {200, 400, 800} x p in {2000, 5000, 8000}
#   beta in {2, 0.9, 0.55} for both
#
# Usage (from the repository root):
#   bash simulation/run_grid.sh [toeplitz|banded] [lowdim|highdim|all]   (defaults: toeplitz all)
#
# If sbatch is available, every cell is submitted as one job via simulation/submit.sh;
# otherwise the cells are run one after another with Rscript.
# Outputs: results/lowdim/, results/highdim/, logs/<regime>_<n>_<p>_<beta>_<setup>.Rout
#---------------------------------------------------------------------------------

if [ ! -f functions/bootstrap_sample.R ]; then
  echo "ERROR: Run this script from the root of the GA-Experiments repository."
  exit 1
fi

SETUP=${1:-toeplitz}
WHICH=${2:-all}

if [ "$SETUP" != "toeplitz" ] && [ "$SETUP" != "banded" ]; then
  echo "Usage: bash simulation/run_grid.sh [toeplitz|banded] [lowdim|highdim|all]"
  exit 1
fi
if [ "$WHICH" != "lowdim" ] && [ "$WHICH" != "highdim" ] && [ "$WHICH" != "all" ]; then
  echo "Usage: bash simulation/run_grid.sh [toeplitz|banded] [lowdim|highdim|all]"
  exit 1
fi

mkdir -p logs

if command -v sbatch > /dev/null 2>&1; then
  USE_SLURM=1
else
  USE_SLURM=0
  echo "NOTE: sbatch not found; running the cells sequentially with Rscript (this takes a long time)."
fi

# Extra sbatch options for the high-dimensional cells, which run in parallel
# (the published runs used 28 CPUs and up to 3 days per cell); adjust to your cluster.
HIGHDIM_SBATCH_OPTS="--cpus-per-task=28 --mem=200G --time=3-00:00:00"

# run_cell <regime> <n> <p> <beta>
run_cell() {
  if [ "$USE_SLURM" -eq 1 ]; then
    if [ "$1" = "highdim" ]; then OPTS=$HIGHDIM_SBATCH_OPTS; else OPTS=""; fi
    sbatch $OPTS --job-name="$1_$2_$3_$4_${SETUP}" simulation/submit.sh "$1" "$2" "$3" "$4" "$SETUP"
  else
    echo "Running $1: n=$2, p=$3, beta=$4, setup=$SETUP"
    Rscript "simulation/$1_sim.R" "$2" "$3" "$4" "$SETUP" > "logs/$1_$2_$3_$4_${SETUP}.Rout" 2>&1 ||
      echo "  FAILED, see logs/$1_$2_$3_$4_${SETUP}.Rout"
  fi
}

for BETA in 2 0.9 0.55; do
  if [ "$WHICH" != "highdim" ]; then
    for N in 200 500 1000 2000 4000 10000; do
      for P in 2 10 30 100; do
        run_cell lowdim "$N" "$P" "$BETA"
      done
    done
  fi
  if [ "$WHICH" != "lowdim" ]; then
    for N in 200 400 800; do
      for P in 2000 5000 8000; do
        run_cell highdim "$N" "$P" "$BETA"
      done
    done
  fi
done
