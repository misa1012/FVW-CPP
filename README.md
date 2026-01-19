
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
Use the provided Python scripts to quantify instability and validate results.
```bash
python tools/python/analyze_instability.py
```
> 🐍 [**Analysis Tools Guide**](docs/analysis_tools.md) - How to use the scripts in `tools/python/`.

---

## 📚 Documentation

The documentation has been reorganized to be more developer-friendly:

- 📘 [**User Guide**](docs/configuration.md): Configuration options, defining perturbation cases.
- 🛠️ [**Developer Guide**](docs/developer_guide.md): Code architecture, adding new C++ tools, build system details.
- 🐍 [**Analysis Tools**](docs/analysis_tools.md): Using the python scripts for post-processing.

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
编辑 `config.json` 设置风速和风机参数。
> 📖 [**配置指南**](docs/configuration.md) - 所有设置项的详细参考。

### 3. 运行
```bash
./fvw_cpp
```

### 4. 分析
使用提供的 Python 脚本量化不稳定性并验证结果。
```bash
python tools/python/analyze_instability.py
```
> 🐍 [**分析工具指南**](docs/analysis_tools.md) - 如何使用 `tools/python/` 中的脚本。

---

## 📚 文档

文档结构经过重组，更加开发者友好：

- 📘 [**用户指南**](docs/configuration.md): 配置选项，定义扰动工况。
- 🛠️ [**开发者指南**](docs/developer_guide.md): 代码架构，添加 C++ 新工具，构建系统细节。
- 🐍 [**分析工具**](docs/analysis_tools.md): 使用 Python 脚本进行后处理。

## 仓库结构

- **`main.cpp`**: 仿真主程序入口。
- **`src/` & `include/`**: C++ 核心库源码。
- **`python/`**: 后处理库 (`fvw`) 和脚本。
- **`data/`**: 风机几何和翼型数据 (如 NREL 5MW)。
- **`tools/`**: 独立的后处理可执行文件 (VTK 网格生成等)。