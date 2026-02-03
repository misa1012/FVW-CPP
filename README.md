
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
Use tutorial configs (strict mode: no default config).
```bash
tutorials/NTNU/config.json
tutorials/NREL/config.json
```
> 📖 [**Configuration Guide**](docs/configuration.md) - Detailed reference for all settings.

### 3. Run
```bash
./fvw_cpp tutorials/NTNU/config.json
```

### 4. Analyze
Use the provided Python scripts to postprocess results.
> 🐍 [**Python Tools**](tools/python/README.md) - How to use the scripts in `tools/python/`.

---

## 📚 Documentation

统一入口：
- 📚 [**Docs Index**](docs/README.md)

## Repository Structure

- **`main.cpp`**: Main simulation entry point.
- **`src/` & `include/`**: C++ Core library source code.
- **`python/`**: Post-processing library (`fvw`) and scripts.
- **`data/`**: Turbine geometry and airfoil data (e.g., NREL 5MW).
- **`tools/`**: Standalone post-processing executables (VTK grid generation, etc.).

---

# FVW-CPP (中文版)

一个风力机自由涡尾迹气动仿真工具，使用 C++17 编写，并配有 Python 后处理工具。

## 🚀 快速开始

### 1. 编译
```bash
mkdir build && cd build
cmake ..
make -j
```

### 2. 配置
使用教程配置（严格模式：不再默认读取根目录配置）。
```bash
tutorials/NTNU/config.json
tutorials/NREL/config.json
```
> 📖 [**配置指南**](docs/configuration.md) - 所有设置项的详细参考。

### 3. 运行
```bash
./fvw_cpp tutorials/NTNU/config.json
```

### 4. 分析
使用提供的 Python 脚本进行后处理。
> 🐍 [**Python 工具**](tools/python/README.md) - 如何使用 `tools/python/` 中的脚本。

---

## 📚 文档

统一入口：
- 📚 [**文档目录**](docs/README.md)

## 仓库结构

- **`main.cpp`**: 仿真主程序入口。
- **`src/` & `include/`**: C++ 核心库源码。
- **`python/`**: 后处理库 (`fvw`) 和脚本。
- **`data/`**: 风机几何和翼型数据 (如 NREL 5MW)。
- **`tools/`**: 独立的后处理可执行文件 (VTK 网格生成等)。
