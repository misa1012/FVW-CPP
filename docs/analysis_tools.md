
# Analysis Tools Guide

The `tools/python/` directory contains standalone scripts for analyzing simulation results. These scripts are designed to work directly with the HDF5 output (`wake.h5`) or CSV probe data.

## Available Tools

### 1. Instability Analysis (`analyze_instability.py`)

**Purpose**: Quantify the onset of tip vortex instabilities (Crow instability, Pairing).

*   **Input**: `results/<CaseName>/wake.h5`
*   **Output**: `results/<CaseName>/post_processing/instability_analysis/`
    *   `radial_instability.png`: Trajectory of tip vortices ($r/x$).
    *   `pairing_distance.png`: Separation distance between adjacent vortices.

**Usage**:
```bash
python tools/python/analyze_instability.py
# (Note: Edit the script directly to change the input filename if needed, default is results/NTNU_Baseline/wake.h5)
```

### 2. Validation with ALM-LES (`compare_fvw_alm.py`)

**Purpose**: Compare mean velocity profiles with high-fidelity ALM-LES data.

*   **Input**: 
    *   FVW: `results/<CaseName>/probe_output.csv`
    *   ALM: Directory path (Simulated or Experimental data)
*   **Output**: `results/<CaseName>/post_processing/comparison_plots/`
    *   Horizontal ($y/R$) and Vertical ($z/R$) velocity profiles at 1D, 2D... 8D.

**Usage**:
```bash
python tools/python/compare_fvw_alm.py
```

### 3. Wake Metrics (`analyze_wake_metrics.py`)

**Purpose**: Compute key wake characteristics like velocity deficit and wake width.

*   **Input**: `results/<CaseName>/probe_output.csv`
*   **Output**: Console output (Metrics).

**Usage**:
```bash
python tools/python/analyze_wake_metrics.py
```

---

# 分析工具指南 (Analysis Tools Guide)

`tools/python/` 目录包含用于分析仿真结果的独立脚本。这些脚本设计为直接处理 HDF5 输出 (`wake.h5`) 或 CSV 探测数据。

## 可用工具

### 1. 不稳定性分析 (`analyze_instability.py`)

**用途**: 量化叶尖涡不稳定性（Crow 不稳定性，配对/蛙跳）的起始点。

*   **输入**: `results/<CaseName>/wake.h5`
*   **输出**: `results/<CaseName>/post_processing/instability_analysis/`
    *   `radial_instability.png`: 叶尖涡的径向轨迹 ($r/x$)。
    *   `pairing_distance.png`: 相邻涡线之间的分离距离。

**用法**:
```bash
python tools/python/analyze_instability.py
# (注意: 如需更改输入文件名，请直接编辑脚本，默认为 results/NTNU_Baseline/wake.h5)
```

### 2. ALM-LES 验证 (`compare_fvw_alm.py`)

**用途**: 将平均速度剖面与高保真 ALM-LES 数据进行对比。

*   **输入**: 
    *   FVW: `results/<CaseName>/probe_output.csv`
    *   ALM: 目录路径 (仿真或实验数据)
*   **输出**: `results/<CaseName>/post_processing/comparison_plots/`
    *   1D, 2D... 8D 处的水平 ($y/R$) 和垂直 ($z/R$) 速度剖面。

**用法**:
```bash
python tools/python/compare_fvw_alm.py
```

### 3. 尾迹指标 (`analyze_wake_metrics.py`)

**用途**: 计算尾迹的关键特征，如速度亏损和尾迹宽度。

*   **输入**: `results/<CaseName>/probe_output.csv`
*   **输出**: 控制台输出 (各项指标)。

**用法**:
```bash
python tools/python/analyze_wake_metrics.py
```
