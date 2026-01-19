
# Configuration Reference (config.json)

FVW-CPP uses `config.json` to control simulation parameters, turbine models, and perturbation cases. This file is the single source of truth for your simulation settings.

## 1. Turbine Configuration

The `turbine` section defines the physical properties of the wind turbine and the inflow conditions.

```json
"turbine": {
    "model": "NREL_5MW",  // Turbine model name. Automatically loads geometry from data/NREL_5MW/
    "windSpeed": 11.4,    // Free stream wind speed (m/s)
    "nSegments": 18,      // Number of radial spanwise segments per blade
    "tsr": 7.0,           // Tip-Speed Ratio (lambda = omega * R / Uinf)
    "rho": 1.225          // Air density (kg/m^3), default is 1.225 if omitted
    
    // Optional Overrides:
    // If "model" is set, these are loaded from data/<model>/turbine_params.json.
    // You can override them here if needed:
    // "rTip": 63.0,      // Rotor radius (m)
    // "rHub": 1.5,       // Hub radius (m)
    // "nBlades": 3,      // Number of blades
    // "hubHeight": 90.0  // Hub height (m)
}
```

### Adding New Models
To add a new turbine model (e.g., `NTNU`):
1. Create `data/NTNU/turbine_params.json`.
2. Add blade geometry to `data/NTNU/blade_geometry.csv`.
3. Set `"model": "NTNU"` in your config.

## 2. Simulation Settings

The `simulation` section controls the time-stepping and numerical parameters of the Free Vortex Wake method.

```json
"simulation": {
    "dt": 0.05,               // Time step size (s)
    "totalTime": 10.0,        // Total simulation duration (s)
    "outputFrequency": 1,     // Write H5 output every N steps (1 = every step)
    "cutoffParam": 0.1,       // Vortex core cut-off parameter (regularization)
    "coreType": "ChordBasedCore", // Core model: "VanGarrel" or "ChordBasedCore"
    "vortexModel": "Constant"     // Vortex decay: "Constant" or "GammaDecay"
}
```

## 3. Perturbation Cases

You can define multiple simulation cases in the `perturbations` array. The main program runs them sequentially.

```json
"perturbations": [
    {
        "name": "case_baseline",   // Output folder name (results/case_baseline/)
        "type": "None",            // Perturbation type
        "amplitude": 0.0,
        "freqFactor": 0.0
    },
    {
        "name": "case_dynamic_pitch",
        "type": "CollectivePitch", 
        "amplitude": 1.0,          // Amplitude in degrees
        "freqFactor": 1.5          // Frequency relative to rotor speed (f = 1.5 * f_rotor)
    }
]
```

### Supported Perturbation Types
| Type | Description | Parameters |
|------|-------------|------------|
| `None` | Baseline constant operation | None |
| `CollectivePitch` | Dynamic sinusoidal collective pitch | `amplitude` (deg), `freqFactor` |
| `AsymmetricStaticPitch` | Static offset applied to one blade | `amplitude` (deg) |

---

# 配置参考手册 (Configuration Reference)

FVW-CPP 使用 `config.json` 来控制仿真参数、风机模型和扰动工况。该文件是所有仿真设置的核心。

## 1. 风机配置 (Turbine Configuration)

`turbine` 部分定义了风机的物理属性和来流条件。

```json
"turbine": {
    "model": "NREL_5MW",  // 风机模型名称。自动从 data/NREL_5MW/ 加载几何信息
    "windSpeed": 11.4,    // 自由来流风速 (m/s)
    "nSegments": 18,      // 每个叶片的径向分段数
    "tsr": 7.0,           // 叶尖速比 (lambda = omega * R / Uinf)
    "rho": 1.225          // 空气密度 (kg/m^3)，如果省略则默认为 1.225
    
    // 可选覆盖项 (Optional Overrides):
    // 如果设置了 "model"，以下参数会从 data/<model>/turbine_params.json 读取。
    // 如果需要，可以在这里强制覆盖：
    // "rTip": 63.0,      // 风轮半径 (m)
    // "rHub": 1.5,       // 轮毂半径 (m)
    // "nBlades": 3,      // 叶片数
    // "hubHeight": 90.0  // 轮毂高度 (m)
}
```

### 添加新模型
要添加一个新的风机模型（例如 `NTNU`）：
1. 创建 `data/NTNU/turbine_params.json`。
2. 添加叶片几何数据到 `data/NTNU/blade_geometry.csv`。
3. 在配置中设置 `"model": "NTNU"`。

## 2. 仿真设置 (Simulation Settings)

`simulation` 部分控制时间步进和自由涡方法的数值参数。

```json
"simulation": {
    "dt": 0.05,               // 时间步长 (s)
    "totalTime": 10.0,        // 仿真总时长 (s)
    "outputFrequency": 1,     // 每 N 步写入一次 H5 输出 (1 = 每步都写)
    "cutoffParam": 0.1,       // 涡核截断参数 (正则化半径)
    "coreType": "ChordBasedCore", // 涡核模型: "VanGarrel" (基于网格长度) 或 "ChordBasedCore" (基于弦长)
    "vortexModel": "Constant"     // 涡衰减模型: "Constant" (恒定) 或 "GammaDecay" (物理衰减)
}
```

## 3. 扰动工况 (Perturbation Cases)

您可以在 `perturbations` 数组中定义多个独立的计算工况。程序会按顺序依次运行它们。

```json
"perturbations": [
    {
        "name": "case_baseline",   // 输出文件夹名称 (生成 results/case_baseline/)
        "type": "None",            // 扰动类型
        "amplitude": 0.0,
        "freqFactor": 0.0
    },
    {
        "name": "case_dynamic_pitch",
        "type": "CollectivePitch", 
        "amplitude": 1.0,          // 幅值 (度)
        "freqFactor": 1.5          // 频率因子 (f = 1.5 * f_rotor)
    }
]
```

### 支持的扰动类型
| 类型 | 描述 | 参数 |
|------|-------------|------------|
| `None` | 基准恒定运行 | None |
| `CollectivePitch` | 动态正弦总距变桨 | `amplitude` (度), `freqFactor` |
| `AsymmetricStaticPitch` | 对单个叶片施加静态变桨偏差 | `amplitude` (度) |
