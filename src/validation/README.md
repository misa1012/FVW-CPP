# Validation Module

此模块包含用于验证 FVW-CPP 模拟结果准确性的工具代码。

## 结构
- **src/**: 包含验证逻辑的实现 (`validate.cpp`, `validate_geometry.cpp`, `validate_position.cpp`).
- **include/**: 包含对应的头文件。

## 功能
主要提供以下验证功能：
1. **Geometry Validation**: 验证叶片几何计算是否符合预期。
2. **Position Validation**: 验证叶片节点在旋转过程中的坐标变换 (DCM 算法) 是否正确。

## 使用方法
目前验证逻辑作为库函数提供。若要在主程序中启用验证，需在 `main` 函数中调用 `fvw::validate` 并传入相应的参数。

示例代码 (在 main.cpp 中):
```cpp
#include "validate.h"

// ... inside main ...
fvw::ValidationParams params;
params.section = "geometry"; // 或 "position"
fvw::validate(params, simParams, turbineParams, geom, pos);
```

未来计划将其封装为独立的测试可执行文件。
