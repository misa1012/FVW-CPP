#!/bin/bash
#SBATCH --job-name=fvw_nrel_15rev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --partition=short

# =============================================================================
# FVW OpenMP Job Submission Script
# =============================================================================

echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $HOSTNAME"
echo "Date:   $(date)"
echo "=================================================="

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
# Navigate to the directory where sbatch was called
cd $SLURM_SUBMIT_DIR

EXECUTABLE="./build/fvw_cpp"
CONFIG="tutorials/NREL/config.json"

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE not found!"
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
