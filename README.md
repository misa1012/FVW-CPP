# FVW-CPP

A free-vortex wake aerodynamic simulation tool for wind turbines, written in C++17 with Python post-processing.

## 🚀 Quick Start

### 1. Build
```bash
mkdir build && cd build
cmake ..
make -j
```

### 2. Configure
Edit `config.json` to set wind speed and turbine parameters.
> 📖 [**Configuration Guide**](docs/configuration.md) - Detailed reference for all settings.

### 3. Run
```bash
./fvw_cpp
```

### 4. Analyze
Use the provided Python package to calculate Cp, Ct, and visualize results.
```bash
cd ../python
python example.py ../build/results/default_baseline/wake.h5
```
> 🐍 [**Python API Guide**](docs/python_api.md) - How to use the `fvw` package.

---

## 📚 Documentation

The documentation has been reorganized to be more developer-friendly:

- 📘 [**User Guide**](docs/configuration.md): Configuration options, defining perturbation cases.
- 🛠️ [**Developer Guide**](docs/developer_guide.md): Code architecture, adding new C++ tools, build system details.
- 🐍 [**Python API**](docs/python_api.md): Using the `fvw` python library for custom analysis.

## Repository Structure

- **`main.cpp`**: Main simulation entry point.
- **`src/` & `include/`**: C++ Core library source code.
- **`python/`**: Post-processing library (`fvw`) and scripts.
- **`data/`**: Turbine geometry and airfoil data (e.g., NREL 5MW).
- **`tools/`**: Standalone post-processing executables (VTK grid generation, etc.).