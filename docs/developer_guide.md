
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

---

# 开发者指南 (Developer Guide)

本指南面向为 FVW-CPP C++ 核心代码做贡献的开发者。

## 架构概览

项目采用模块化架构，将物理模型、数据处理和仿真控制分离开来。

### 目录结构
- **`src/`**: 共享库代码 (`libfvw_core.a`)。
    - **`core/`**: 基础数据结构 (`Wake`, `Blade`, `Node`)。
    - **`models/`**: 气动模型 (`_geometry`, `Airfoils`, `induced_velocity`)。
    - **`io/`**: 输入/输出处理 (`ConfigLoader`, `H5Writer`)。
    - **`simulation/`**: 高层驱动器 (`SimulationRunner`)。
- **`include/`**: 公共头文件 (与 `src/` 结构对应)。
    - **`core/`**: `geometry.h`, `wake.h` 等。
    - **`models/`**: `airfoil.h`, `performance.h`, `bem.h`。
    - **`io/`**: `config.h`, `postprocess.h`。
    - **`simulation/`**: `simulation_runner.h`。
- **`tools/`**: 辅助后处理可执行文件 (`postprocess_grid` 等)。
- **`main.cpp`**: 仿真程序主入口。

### 关键类

#### `SimulationRunner` (`src/simulation/simulation_runner.cpp`)
单个仿真算例的指挥官。
- **角色**: 初始化模型、运行时间循环、处理 I/O。
- **生命周期**: `initialize()` -> `run()` -> `finalize()`。

#### `Wake` (`src/core/wake.cpp`)
管理拉格朗日尾迹的表示。
- 存储涡线和粒子。
- 处理尾迹的对流 (Convection) 和脱落 (Shedding) 步骤。

#### `ConfigLoader` (`src/io/config.cpp`)
解析 `config.json` 并加载模型特定的默认参数 (例如 `data/NREL_5MW/turbine_params.json`)。

## 构建系统

项目使用 CMake (3.10+)。

### 编译目标 (Targets)
- **`fvw_core`**: 包含所有共享逻辑的静态库。
- **`fvw_cpp`**: 链接到 `fvw_core` 的主可执行文件。
- **`postprocess_XXX`**: 链接到 `fvw_core` 的小工具。

### 构建配置选项

您可以使用 CMake 选项自定义构建过程：

| 选项 | 默认值 | 描述 |
|--------|---------|-------------|
| `BUILD_FVW_CPP` | `ON` | 主仿真可执行文件。 |
| `BUILD_TOOL_GRID` | `ON` | VTK 网格生成后处理工具。 |
| `BUILD_TOOL_PROBES` | `OFF` | 点探测器采样工具。 |
| `BUILD_TOOL_DOWNSTREAM` | `OFF` | 下游速度剖面计算器。 |
| `BUILD_VALIDATION` | `OFF` | 编译验证库 (检查几何/位置数学)。 |
| `ENABLE_PROFILING` | `OFF` | 启用 `gprof` 性能分析标志。 |

**示例**: 仅编译主仿真程序，不编译任何工具：
```bash
cmake .. -DBUILD_TOOL_GRID=OFF
make
```
