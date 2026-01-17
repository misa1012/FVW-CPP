# Python Post-Processing API

The `python/` directory contains a modular Python package `fvw` for analyzing simulation results stored in HDF5 format.

## Overview

Instead of monolithic scripts, we provide a reusable library that handles data loading and aerodynamic calculations.

- **`fvw.io`**: Handles HDF5 reading, configuration parsing, and geometry fallbacks.
- **`fvw.calc`**: Implements Blade Element Theory (BET) integration to compute Thrust, Torque, Power, $C_P$, and $C_T$.

## Usage Guide

### 1. Setup
Ensure you have the required dependencies:
```bash
pip install h5py numpy matplotlib
```

### 2. Loading Data
The `WakeReader` class standardizes data access.

```python
from fvw.io import WakeReader

# Initialize with path to wake.h5
reader = WakeReader("results/case_baseline/wake.h5")

# Get simulation parameters
config = reader.read_config()
print(f"Rotor Speed: {config['omega']} rad/s")

# Get available timesteps
timesteps = reader.list_timesteps()
```

### 3. Computing Performance
The `PerformanceCalculator` computes integrated loads.

```python
from fvw.calc import PerformanceCalculator

calc = PerformanceCalculator(reader)

# Compute time series, skipping initial transients (start_rev=5)
# This calculates Thrust/Power at every timestep by integrating blade loads
res = calc.compute_time_series(start_rev=5.0)

# Access results
time = res['time']
cp = res['CP']
ct = res['CT']

print(f"Mean Cp: {cp.mean():.4f}")
```

## Analysis Tool
See `python/process_aero.py` for a complete runnable analysis tool that plots $C_P$, $C_T$, Power, and Thrust over time.

```bash
cd python
python process_aero.py ../results/case_baseline/wake.h5
```

### Spanwise Analysis
Use `python/process_spanwise.py` to inspect coefficients ($C_l, C_d, \alpha$) along the blade span for a specific timestep.

```bash
# Analyze the last valid timestep (default)
python process_spanwise.py ../results/case_baseline/wake.h5

# Analyze a specific timestep (e.g., 100)
python process_spanwise.py ../results/case_baseline/wake.h5 --step 100
```
