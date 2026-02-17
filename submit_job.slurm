#!/bin/bash
#SBATCH --job-name=fvw_ntnu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=08:00:00
#SBATCH --output=/dev/null
#SBATCH --partition=short

# =============================================================================
# FVW OpenMP Job Submission Script
# =============================================================================

# --- 1. Environment Setup ---
# Uncomment and modify specific module loads for your HPC system
module purge
module load GCC/11.2.0
module load HDF5/1.14.0-GCC-11.2.0-serial
module load CMake/3.22.1-GCCcore-11.2.0

# --- 2. OpenMP Configuration ---
# Crucial for performance on high-core-count nodes
# Set number of threads to match the requested CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Thread Binding: Prevents OS from migrating threads between cores
export OMP_PROC_BIND=close 
export OMP_PLACES=cores

# Debug: Print OpenMP environment to log
export OMP_DISPLAY_ENV=true 

# --- 3. Execution ---
# Resolve executable robustly. On this cluster, Slurm may execute a spooled
# copy of the script, so BASH_SOURCE can point to /var/spool/... and break
# relative project-root logic.
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"

EXEC_CANDIDATES=(
    "${FVW_EXECUTABLE:-}"
    "/data/marine-cfd/shug8104/FVW-CPP/build/fvw_cpp"
    "/home/shug8104/data/FVW-CPP/build/fvw_cpp"
    "$SLURM_SUBMIT_DIR/build/fvw_cpp"
    "$SCRIPT_DIR/build/fvw_cpp"
)

EXECUTABLE=""
for cand in "${EXEC_CANDIDATES[@]}"; do
    [[ -n "$cand" ]] || continue
    if [[ -x "$cand" ]]; then
        EXECUTABLE="$cand"
        break
    fi
done
CONFIG="${CONFIG:-tutorials/NTNU/config.json}"

# Route stdout/stderr to the results dir for the case (same place as simulation.log)
OUTDIR=$(python - <<'PY'
import json, os
cfg = os.environ.get("CONFIG", "tutorials/NTNU/config.json")
with open(cfg) as f:
    data = json.load(f)
case = data.get("caseName")
root = data.get("outputRoot")
model = data.get("turbine", {}).get("model", "case")
if root and case:
    print(os.path.join(root, case))
elif root:
    print(os.path.join(root, model))
elif case:
    print(os.path.join("results", case))
else:
    print(os.path.join("results", model))
PY
)
case "$OUTDIR" in
    /*) ;; # absolute path, keep as-is
    *) OUTDIR="$SLURM_SUBMIT_DIR/$OUTDIR" ;;
esac
mkdir -p "$OUTDIR"
OUTFILE="$OUTDIR/slurm-$SLURM_JOB_ID.out"
exec > "$OUTFILE" 2>&1

echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $HOSTNAME"
echo "Date:   $(date)"
echo "=================================================="

if [ -z "$EXECUTABLE" ] || [ ! -x "$EXECUTABLE" ]; then
    echo "Error: fvw_cpp executable not found in candidates:"
    for cand in "${EXEC_CANDIDATES[@]}"; do
        [[ -n "$cand" ]] && echo "  $cand"
    done
    echo "Please compile before submitting:"
    echo "  mkdir -p build && cd build"
    echo "  cmake .. -DCMAKE_BUILD_TYPE=Release"
    echo "  make -j$SLURM_CPUS_PER_TASK"
    exit 1
fi

echo "Running Simulation with $OMP_NUM_THREADS threads..."
time $EXECUTABLE $CONFIG

echo "=================================================="
echo "Job Completed at $(date)"
echo "=================================================="
