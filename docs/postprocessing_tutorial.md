# Post-Processing Tutorial (NREL & NTNU)

This project uses two post-processing tools:

- `postprocess_grid`: export VTK (2D slice or 3D volume)
- `postprocess_sampling`: sample induced velocity at specified points/lines, instant or average

Note: The correct turbine name is **NREL**, not "NERL".

## Where to put sampling.json

Recommended:
- Reusable templates: `docs/tutorials/<CASE>/sampling.json`
- For traceability: copy the actual file to `results/<case_name>/sampling.json`

## Build (post-processing tools)

```bash
mkdir -p build
cd build
cmake .. -DBUILD_TOOL_GRID=ON -DBUILD_TOOL_SAMPLING=ON
make -j
```

## Tutorial A: NREL_5MW

Example files:
- `tutorials/NREL/config.json`
- `tutorials/NREL/sampling.json`

Run simulation:
```bash
./build/fvw_cpp tutorials/NREL/config.json
```

Sampling (average over a time window):
```bash
./build/postprocess_sampling tutorials/NREL/sampling.json
```

VTK export (3D or slice):
```bash
./build/postprocess_grid results/NREL_5MW/wake.h5
```

Optional VTK parameters:
```
postprocess_grid <wake.h5> [out.vtk] [config.json] [timestep] [grid_res] [slice_mode] [x_start_factor] [x_end_factor]
```
- `slice_mode`: 0 = full 3D, 1 = horizontal plane (Z fixed), 2 = vertical plane (Y fixed)
- `grid_res`: grid spacing in meters
- `x_start_factor`, `x_end_factor`: bounds in rotor diameters (D)

## Tutorial B: NTNU

Example files:
- `tutorials/NTNU/config.json`
- `tutorials/NTNU/sampling.json`

Run simulation:
```bash
./build/fvw_cpp tutorials/NTNU/config.json
```

Sampling:
```bash
./build/postprocess_sampling tutorials/NTNU/sampling.json
```

VTK export:
```bash
./build/postprocess_grid results/NTNU_0_25/wake.h5
```

## sampling.json format

Minimal example:
```json
{
  "wake_h5": "results/<case>/wake.h5",
  "config": "tutorials/NTNU/config.json",
  "output_csv": "results/<case>/sampling_output.csv",
  "mode": "instant",
  "instant": { "step": 100 },
  "points": [
    { "name": "p0", "x": 1.0, "y": 0.0, "z": 0.0 }
  ]
}
```

Average example:
```json
{
  "wake_h5": "results/<case>/wake.h5",
  "config": "tutorials/NTNU/config.json",
  "output_csv": "results/<case>/sampling_output.csv",
  "mode": "average",
  "units": "m",
  "average": { "start_step": 100, "end_step": 200, "stride": 5 },
  "lines": [
    { "name": "centerline", "start": [1.0, 0.0, 90.0], "end": [5.0, 0.0, 90.0], "n": 20 }
  ]
}
```

Notes:
- `units`: `"m"` (default) or `"D"` (rotor diameter)
- `points`: list of discrete points
- `lines`: line sampling with `n` points (inclusive of endpoints)

---

# 后处理教程（NREL & NTNU）

本项目后处理主要包含两类工具：
- `postprocess_grid`：导出 VTK（2D 切片或 3D 体）
- `postprocess_sampling`：按指定位置采样，支持瞬时或平均

## 建议放置 sampling.json
- 可复用模板：`tutorials/<CASE>/sampling.json`
- 结果可追溯：复制一份到 `results/<case_name>/sampling.json`

## 编译后处理工具
```bash
mkdir -p build
cd build
cmake .. -DBUILD_TOOL_GRID=ON -DBUILD_TOOL_SAMPLING=ON
make -j
```

## 教程 A：NREL_5MW

示例文件：
- `tutorials/NREL/config.json`
- `tutorials/NREL/sampling.json`

运行：
```bash
./build/fvw_cpp tutorials/NREL/config.json
./build/postprocess_sampling tutorials/NREL/sampling.json
./build/postprocess_grid results/NREL_5MW/wake.h5
```

## 教程 B：NTNU

示例文件：
- `tutorials/NTNU/config.json`
- `tutorials/NTNU/sampling.json`

运行：
```bash
./build/fvw_cpp tutorials/NTNU/config.json
./build/postprocess_sampling tutorials/NTNU/sampling.json
./build/postprocess_grid results/NTNU_0_25/wake.h5
```

## sampling.json 格式

**瞬时采样：**
```json
{
  "wake_h5": "results/<case>/wake.h5",
  "config": "tutorials/NTNU/config.json",
  "output_csv": "results/<case>/sampling_output.csv",
  "mode": "instant",
  "instant": { "step": 100 },
  "points": [
    { "name": "p0", "x": 1.0, "y": 0.0, "z": 0.0 }
  ]
}
```

**时间平均：**
```json
{
  "wake_h5": "results/<case>/wake.h5",
  "config": "tutorials/NTNU/config.json",
  "output_csv": "results/<case>/sampling_output.csv",
  "mode": "average",
  "units": "m",
  "average": { "start_step": 100, "end_step": 200, "stride": 5 },
  "lines": [
    { "name": "centerline", "start": [1.0, 0.0, 90.0], "end": [5.0, 0.0, 90.0], "n": 20 }
  ]
}
```

备注：
- `units`：`"m"`（默认）或 `"D"`（以转子直径归一化）
- `points`：离散采样点
- `lines`：线采样（自动生成 n 个点）
