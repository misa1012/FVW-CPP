
# Research Plan: Tip Vortex Instability Mechanisms

## 1. Objective
Investigate the onset and evolution of tip vortex instabilities (specifically **Crow Instability** and **Pairing/Leapfrogging**) using the FVW-CPP solver.
The goal is to quantify how parameters like **Tip Speed Ratio (TSR)**, **Vortex Core Size (Cutoff)**, and **Perturbation Amplitude** affect the instability growth rate and breakdown location.

## 2. Methodology
The study utilizes the Inviscid Free Vortex Wake (FVW) method. Since the model is naturally stable, **external perturbations** must be introduced to trigger physically meaningful instabilities.

### Key Variables
*   **Cutoff Parameter (`cutoffParam`)**: Controls the vortex core radius.
    *   *Range*: 0.05 (Fine) to 0.25 (Coarse).
    *   *Hypothesis*: Smaller cutoff $\to$ Faster short-wave growth (Pairing).
*   **Perturbation**:
    *   *Type*: `AsymmetricStaticPitch` (Geometry) or `CollectivePitch` (Transient).
    *   *Amplitude*: $0.5^\circ \sim 2.0^\circ$.
*   **TSR**: Tip Speed Ratio (Affects helix pitch and mutual inductance).

## 3. Workflow & Tools
All analysis scripts are located in `tools/python/`.

### Step 1: Simulation Setup
使用教程配置（示例）：`tutorials/NTNU/config.json`
```json
"perturbations": [
    {
        "type": "AsymmetricStaticPitch", 
        "amplitude": 1.0, 
        "freqFactor": 0.0
    }
],
"simulation": {
    "cutoffParam": 0.1
}
```

### Step 2: Execution
Run the simulation on the cluster:
```bash
sbatch submit_job.slurm
```
(Ensure `outputFrequency` is small enough, e.g., 10, to capture dynamics).

### Step 3: Post-Processing & Analysis

#### A. Instability Analysis (Primary Tool)
Quantify radial wobble and pairing.
*   **Script**: `tools/python/analyze_instability.py`
*   **Usage**:
    ```bash
    python tools/python/analyze_instability.py
    ```
*   **Outputs** (`results/NTNU_Baseline/post_processing/instability_analysis/`):
    *   `radial_instability.png`: Visualizes $r/x$ trajectory. Look for oscillating growth (Crow) or straight lines (Stable).
    *   `pairing_distance.png`: Distance between adjacent vortices. Look for sudden drops (Leapfrogging).

#### B. Validation with ALM-LES
Compare mean velocity profiles with high-fidelity LES data.
*   **Script**: `tools/python/compare_fvw_alm.py`
*   **Usage**:
    ```bash
    python tools/python/compare_fvw_alm.py
    ```
*   **Requires**: `results/NTNU_Baseline/probe_output.csv` (Use postprocessing sampling tool to generate).
*   **Outputs**: `results/NTNU_Baseline/post_processing/comparison_plots/`

#### C. Wake Metrics
Compute deficit and expansion trends.
*   **Script**: `tools/python/analyze_wake_metrics.py` (Optional)

## 4. File Cleanup Strategy
To save space and keep the workspace clean:
*   **Keep**: `wake.h5` (Core data), `config_*.json` (Record of settings), `analysis_*` (Plots).
*   **Delete**: 
    *   `probe_output_partial.csv` (Intermediate debug file).
    *   `slurm-*.out` (Old logs, unless containing errors).
    *   `debug_*.txt` (If any).

## 5. Next Steps Experiment Matrix

| Experiment ID | Cutoff | Perturbation | Objective |
| :--- | :--- | :--- | :--- |
| **Exp-01 (Baseline)** | 0.1 | None | Confirm absolute stability (Straight lines). |
| **Exp-02 (Low Noise)**| 0.25 | None | Verify removal of numerical noise/ejections. |
| **Exp-03 (Forced)**   | 0.25 | 1.0 deg Pitch | Trigger & Measure Crow Instability Growth Rate. |
| **Exp-04 (Fine Core)**| 0.05 | 1.0 deg Pitch | Investigate Short-wave (Pairing) sensitivity. |

---

# 研究计划：叶尖涡不稳定性机制

## 1. 目标
利用 FVW-CPP 求解器研究叶尖涡不稳定性的起始和演化（特别是 **Crow 不稳定性** 和 **涡配对/蛙跳 (Pairing/Leapfrogging)**）。
目标是量化 **叶尖速比 (TSR)**、**涡核尺寸 (Cutoff)** 和 **扰动幅值** 等参数如何影响不稳定性的增长率和破碎位置。

## 2. 方法论
本研究使用无粘自由涡尾迹 (FVW) 方法。由于该模型本身是自然稳定的，必须引入 **外部扰动** 才能触发具有物理意义的不稳定性。

### 关键变量
*   **截断参数 (`cutoffParam`)**: 控制涡核半径。
    *   *范围*: 0.05 (精细) 到 0.25 (粗糙)。
    *   *假设*: 较小的 cutoff $\to$ 导致短波增长更快 (Pairing)。
*   **扰动**:
    *   *类型*: `AsymmetricStaticPitch` (非对称静态变桨) 或 `CollectivePitch` (周期性总距变桨)。
    *   *幅值*: $0.5^\circ \sim 2.0^\circ$。
*   **TSR**: 叶尖速比 (影响螺旋螺距和互感强度)。

## 3. 工作流与工具
所有分析脚本均位于 `tools/python/`。

### 步骤 1: 仿真设置
编辑教程配置（示例）：`tutorials/NTNU/config.json`
```json
"perturbations": [
    {
        "type": "AsymmetricStaticPitch", 
        "amplitude": 1.0, 
        "freqFactor": 0.0
    }
],
"simulation": {
    "cutoffParam": 0.1
}
```

### 步骤 2: 执行
在集群上运行仿真：
```bash
sbatch submit_job.slurm
```
(确保 `outputFrequency` 足够小，例如 10，以捕捉动态过程)。

### 步骤 3: 后处理与分析

#### A. 不稳定性分析 (核心工具)
量化径向摆动和配对现象。
*   **脚本**: `tools/python/analyze_instability.py`
*   **用法**:
    ```bash
    python tools/python/analyze_instability.py
    ```
*   **输出** (`results/NTNU_Baseline/post_processing/instability_analysis/`):
    *   `radial_instability.png`: 可视化 $r/x$ 轨迹。寻找震荡增长 (Crow) 或直线 (Stable)。
    *   `pairing_distance.png`: 相邻涡线间的距离。寻找突然下降 (Leapfrogging)。

#### B. ALM-LES 验证
将平均速度剖面与高保真 LES 数据进行对比。
*   **脚本**: `tools/python/compare_fvw_alm.py`
*   **用法**:
    ```bash
    python tools/python/compare_fvw_alm.py
    ```
*   **依赖**: `results/NTNU_Baseline/probe_output.csv` (使用后处理采样工具生成)。
*   **输出** (`results/NTNU_Baseline/post_processing/comparison_plots/`):

#### C. 尾迹指标
计算亏损和膨胀趋势。
*   **脚本**: `tools/python/analyze_wake_metrics.py` (可选)

## 4. 文件清理策略
为了节省空间并保持工作区整洁：
*   **保留**: `wake.h5` (核心数据), `config_*.json` (设置记录), `analysis_*` (图表)。
*   **删除**: 
    *   `probe_output_partial.csv` (中间调试文件)。
    *   `slurm-*.out` (旧日志，除非含报错)。
    *   `debug_*.txt` (如果有)。

## 5. 下一步实验矩阵

| 实验 ID | Cutoff | 扰动 | 目标 |
| :--- | :--- | :--- | :--- |
| **Exp-01 (基准)** | 0.1 | 无 | 确认绝对稳定性 (直线)。 |
| **Exp-02 (低噪)** | 0.25 | 无 | 验证数值噪音/弹射点的去除效果。 |
| **Exp-03 (强迫)** | 0.25 | 1.0度 变桨 | 触发并测量 Crow 不稳定性的增长率。 |
| **Exp-04 (精细核)**| 0.05 | 1.0度 变桨 | 考察短波 (Pairing) 的敏感性。 |
