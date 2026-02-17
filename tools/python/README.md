# FVW-CPP Python Analysis Tools

This directory contains the post-processing and analysis suite for the FVW solver. Tools are organized by **function**, plus a **legacy** area for old or one-off scripts.

## 1. Directory Structure

```
tools/python/
├── aerodynamics/
│   ├── single_case/         # Single-case metrics/plots (Cp/Ct/Power/Thrust/AoA/Cl/Cd/Fn/Ft)
│   └── compare_alm/         # Comparisons with ALM-LES
├── downstream_velocity/     # Wake profiles, deficit, REWS, metrics
├── wake_visualization/      # FVW / ALM / LES visualization
├── studies/                 # Study-specific workflows (e.g., cutoff)
│   └── cutoff/
└── legacy/                  # Old / one-off scripts
    └── spatial_visualization/
```

## 2. Usage Guide

### Prerequisites
Activate the conda environment with necessary libraries (numpy, pandas, matplotlib, pyvista, h5py):
```bash
conda run -n post python3 ...
```

### A. Aerodynamics

#### 1) Single Case Metrics (plots on demand)
Generate selected plots for a single case:
```bash
conda run -n post python tools/python/aerodynamics/single_case/analyze_case.py \
    results/CASE_NAME/wake.h5 --plots power thrust cp ct
```

Spanwise metrics at last timestep (or time-average):
```bash
conda run -n post python tools/python/aerodynamics/single_case/analyze_case.py \
    results/CASE_NAME/wake.h5 --plots aoa cl cd fn ft --timestep last
```

#### 2) Compare with ALM-LES
Cp/Ct comparison:
```bash
conda run -n post python tools/python/aerodynamics/compare_alm/compare.py \
    results/CASE_NAME/aerodynamic_performance/wake_perf.csv
```

AoA comparison:
```bash
conda run -n post python tools/python/aerodynamics/compare_alm/compare_aoa_alm.py \
    results/CASE_NAME/wake.h5 --alm-dir /path/to/ALM_CASE
```

### B. Downstream Velocity (Wake Deficit)

**Extract Profiles**: Extracts velocity lines at 1D, 2D...8D downstream using high-res sampling.
```bash
conda run -n post python tools/python/downstream_velocity/extract_profiles.py \
    results/CASE_NAME/wake.h5 \
    -o results/CASE_NAME/probe_output.csv
```

**Compare with ALM-LES**: Generates deficit plots ($1 - U/U_{inf}$) comparing FVW and ALM.
```bash
conda run -n post python tools/python/downstream_velocity/compare_alm.py \
    results/CASE_NAME/probe_output.csv \
    -r /path/to/alm/les/case
```

**Calculate Metrics**: Prints summary table of max deficit and error % vs ALM.
```bash
conda run -n post python tools/python/downstream_velocity/calculate_metrics.py \
    results/CASE_NAME/probe_output.csv \
    -r /path/to/alm/les/case
```

**REWS**: Generates YZ-plane deficit maps for FVW with REWS calculation.
```bash
conda run -n post python tools/python/downstream_velocity/calculate_rews.py \
    results/ntnu_fixed_fs10_v2/wake.h5
```

### C. Wake Visualization

**Visualize ALM-LES YZ Wakes**: Generates YZ-plane velocity deficit maps from ALM VTK/OpenFOAM data.
```bash
conda run -n post python tools/python/wake_visualization/visualize_alm_wake.py \
    --alm_dir /path/to/ALM_CASE
```

**Visualize FVW 3D Field**:
```bash
conda run -n post python tools/python/wake_visualization/visualize_fvw_wake.py \
    --results_dir results/CASE_NAME
```
(Note: Ensure `postprocess_grid` is compiled in `build/`.)

### D. Studies (Optional)

**Cutoff Study** (VTK + deficit plots):
```bash
conda run -n post python tools/python/studies/cutoff/visualize.py \
    --study-dir /path/to/cutoff_study \
    --project-root /path/to/FVW-CPP \
    --config-template /path/to/tutorials/NTNU/config.json
```

**Compare Cutoff Study** (requires probe_output.csv in each case):
```bash
conda run -n post python tools/python/studies/cutoff/compare.py \
    --study-dir /path/to/cutoff_study \
    --alm-dir /path/to/alm/case
```

**Cutoff Postprocess (Cp/Ct + spanwise)**:
```bash
# Multi-case
conda run -n post python tools/python/studies/cutoff/postprocess_cutoff.py \
    --study-dir /path/to/cutoff_study --plots cp ct aoa cl fn ft

# Single case
conda run -n post python tools/python/studies/cutoff/postprocess_cutoff.py \
    --case-dir /path/to/CASE_DIR --plots cp ct
```

**Cutoff Wake Plots (velocity deficit / vorticity)**:
```bash
# Multi-case
conda run -n post python tools/python/studies/cutoff/wake_compare_cutoff.py \
    --study-dir /path/to/cutoff_study --ppd 40 --xlim -0.5 2.0 --ylim -0.5 0.5

# Single case
conda run -n post python tools/python/studies/cutoff/wake_compare_cutoff.py \
    --case-dir /path/to/CASE_DIR --ppd 40 --xlim -0.5 2.0 --ylim -0.5 0.5
```

---
**Tips:**
- Use `-h` on any script to see available options.
- Ensure `probe_output.csv` exists before running comparison scripts.
