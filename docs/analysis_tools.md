
# Analysis Tools Guide / 分析工具指南

本页合并了原 `analysis_tools.md` 与 `analysis_tool_guide.md`，给出统一入口与推荐用法。

## 入口 / Entry
- Python 工具总览：`tools/python/README.md`
- 后处理教程：`docs/postprocessing_tutorial.md`

## 1. 可视化工具 / Visualization

### FVW 结果可视化
脚本：`tools/python/wake_visualization/visualize_fvw_wake.py`  
用途：生成尾迹速度亏损、涡量等可视化图  
示例：
```bash
conda run -n post python tools/python/wake_visualization/visualize_fvw_wake.py \
  --results_dir results/CASE_NAME
```

### ALM/LES 结果可视化
脚本：`tools/python/wake_visualization/visualize_alm_wake.py`  
用途：读取 ALM-LES 输出进行对比展示  
示例：
```bash
conda run -n post python tools/python/wake_visualization/visualize_alm_wake.py --alm_dir /path/to/ALM_CASE
```

## 2. 定量分析 / Quantitative

### 气动性能 (Cp/Ct)
目录：`tools/python/aerodynamics/`  
主要脚本：
- `single_case/analyze_case.py`：单 case 计算/绘图（Cp/Ct/Power/Thrust/AoA/Cl/Cd/Fn/Ft）
- `compare_alm/compare.py`：Cp/Ct 与 ALM 参考对比
- `compare_alm/compare_aoa_alm.py`：AoA 与 ALM 参考对比

示例：
```bash
conda run -n post python tools/python/aerodynamics/single_case/analyze_case.py \
  results/CASE/wake.h5 --plots power thrust cp ct

conda run -n post python tools/python/aerodynamics/compare_alm/compare.py \
  results/CASE/aerodynamic_performance/wake_perf.csv
```

### 下游速度剖面
目录：`tools/python/downstream_velocity/`  
用途：提取剖面、对比 ALM、计算指标  

示例：
```bash
conda run -n post python tools/python/downstream_velocity/extract_profiles.py results/CASE/wake.h5
conda run -n post python tools/python/downstream_velocity/compare_alm.py results/CASE/probe_output.csv -r /path/to/ALM
```

### Cutoff Study
目录：`tools/python/studies/cutoff/`
```bash
conda run -n post python tools/python/studies/cutoff/visualize.py \
  --study-dir /path/to/cutoff_study \
  --project-root /path/to/FVW-CPP \
  --config-template tutorials/NTNU/config.json
```

## 3. 不稳定性分析 / Instability
脚本：`tools/python/downstream_velocity/analyze_instability.py`  
用途：诊断尾迹不稳定性  
示例：
```bash
conda run -n post python tools/python/downstream_velocity/analyze_instability.py
```
