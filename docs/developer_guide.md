# Developer Guide

This guide is for developers contributing to the C++ core of FVW-CPP.

## Architecture Overview

The project follows a modular architecture separating physics models, data handling, and simulation control.

### Directory Structure
- **`src/`**: Shared library code (`libfvw_core.a`).
    - **`core/`**: Fundamental data structures (`Wake`, `Blade`, `Node`).
    - **`models/`**: Aerodynamic models (`Geometry`, `Airfoils`, `induced_velocity`).
    - **`io/`**: Input/Output handling (`ConfigLoader`, `H5Writer`).
    - **`simulation/`**: High-level driver (`SimulationRunner`).
- **`include/`**: Public headers (Mirrors `src/` structure).
    - **`core/`**: `geometry.h`, `wake.h`, etc.
    - **`models/`**: `airfoil.h`, `performance.h`, `bem.h`
    - **`io/`**: `config.h`, `postprocess.h`
    - **`simulation/`**: `simulation_runner.h`
- **`tools/`**: Auxiliary post-processing executables (`postprocess_grid`, etc.).
- **`main.cpp`**: The primary simulation entry point.

### Key Classes

#### `SimulationRunner` (`src/simulation/simulation_runner.cpp`)
The orchestrator of a single simulation case.
- **Role**: Initializes models, runs the time loop, and handles I/O.
- **Life Cycle**: `initialize()` -> `run()` -> `finalize()`.

#### `Wake` (`src/core/wake.cpp`)
Manages the Lagrangian wake representation.
- Stores vortex filaments and particles.
- Handles convection and shedding steps.

#### `ConfigLoader` (`src/io/config.cpp`)
Parses `config.json` and loads model-specific defaults (e.g., `data/NREL_5MW/turbine_params.json`).

## Build System

The project uses CMake (3.10+).

### Targets
- **`fvw_core`**: Static library containing all shared logic.
- **`fvw_cpp`**: Main executable linked against `fvw_core`.
- **`postprocess_XXX`**: Tools linked against `fvw_core`.

### Build Configurations

You can customize the build using CMake options:

| Option | Default | Description |
|--------|---------|-------------|
| `BUILD_FVW_CPP` | `ON` | Main simulation executable. |
| `BUILD_TOOL_GRID` | `ON` | VTK grid generation post-processing tool. |
| `BUILD_TOOL_PROBES` | `OFF` | Point probe sampling tool. |
| `BUILD_TOOL_DOWNSTREAM` | `OFF` | Downstream velocity profile calculator. |
| `BUILD_VALIDATION` | `OFF` | Compile validation library (checking geometry/position math). |
| `ENABLE_PROFILING` | `OFF` | Enable `gprof` profiling flags. |

**Example**: Build only the main simulation without any tools:
```bash
cmake .. -DBUILD_TOOL_GRID=OFF
make
```

