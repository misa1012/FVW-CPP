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
- **`include/`**: Public headers.
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

### Adding New Tools
To add a new C++ tool:
1. Create source file in `tools/new_tool.cpp`.
2. Add executable in `CMakeLists.txt`:
   ```cmake
   add_executable(new_tool tools/new_tool.cpp)
   target_link_libraries(new_tool PRIVATE fvw_core)
   ```
