# FVW-CPP 用户手册：原理与实现指南

本文档是一份完整的参考手册，旨在帮助用户从理论与实践两个维度理解 FVW-CPP 项目。

本文档包含三个部分：
1.  **物理原理 (Theory)**：自由尾迹方法的基本假设与数学公式。
2.  **代码实现对照 (Map)**：物理公式在 C++ 代码中的具体实现位置。
3.  **运行指南 (Run Guide)**：如何配置与运行仿真。

---

## Part 1: 物理原理 (Theory)

本部分详细介绍了自由尾迹方法（Free Vortex Wake, FVW）的物理核心。重点解释叶片气动力计算、相对速度合成以及尾迹演化过程。

### 1.1 核心思想

自由尾迹方法是一种势流理论方法，用于模拟风力机或螺旋桨的非定常气动特性。其核心假设是流场是无粘、不可压缩的，通过离散的“涡线”（Vortex Filaments）来模拟叶片产生的升力和尾迹对流场的诱导作用。

**“自由”的含义**：尾迹的几何形状不是预先固定的（如刚性尾迹模型），而是随着流场中的局部速度自由畸变（对流），从而捕捉尾迹膨胀、卷起和相互干扰等非线性现象。

### 1.2 坐标系与速度三角形

在任意时刻 $t$，叶片某一截面（翼型）主要受到三种速度的影响，它们的矢量和构成了**相对风速 ($\vec{V}_{rel}$)**：

$$ \vec{V}_{rel} = (\vec{V}_{\infty} + \vec{V}_{ind}) - \vec{V}_{motion} $$

1.  **$\vec{V}_{\infty}$ (自由来流)**：远前方吹来的风速。
2.  **$\vec{V}_{ind}$ (诱导速度)**：由全场涡系在叶片处产生的诱导速度（通过 Biot-Savart 定律计算）。
3.  **$\vec{V}_{motion}$ (运动速度)**：叶片自身的运动速度。
    *   **物理直观**：如果你骑自行车（叶片动）向东，即使没有风，你也会感觉到风从东边吹来（$-\vec{V}_{motion}$）。

### 1.3 气动力计算 (Kutta-Joukowski)

一旦得到了相对风速 $\vec{V}_{rel}$，我们可以计算攻角 $\alpha$ 并查表得到升力系数 $C_l$。根据 Kutta-Joukowski 定理，升力与环量 $\Gamma$ 相关：

$$ L = \rho |\vec{V}_{rel}| \Gamma \implies \Gamma = \frac{1}{2} |\vec{V}_{rel}| c C_l(\alpha) $$

这构成了一个非线性迭代过程：$\Gamma$ 决定 $\vec{V}_{ind}$，$\vec{V}_{ind}$ 决定 $\alpha$，$\alpha$ 决定 $C_l$，而 $C_l$ 又决定新的 $\Gamma$。

### 1.4 尾迹演化 (Wake Convection)
根据 Helmholtz 涡定理，涡线随流体运动。尾迹节点的更新遵循拉格朗日法：

$$ \vec{x}(t + \Delta t) = \vec{x}(t) + [\vec{V}_{\infty} + \vec{V}_{ind}(\vec{x})] \cdot \Delta t $$

---

## Part 2: 代码实现对照 (Implementation Map)

本部分帮助您在 C++ 代码中找到上述物理原理的具体实现。

### 2.1 核心函数速查表

| 物理过程 | 核心函数 |所在文件 | 说明 |
| :--- | :--- | :--- | :--- |
| **1. 叶片运动** | `computePositions` | `src/core/position.cpp` | 根据转速 $\Omega$ 更新叶片节点坐标 $\vec{r}(t)$ |
| **2. 运动速度** | `computeVelICS` | `src/core/velocity.cpp` | 计算 $\vec{V}_{motion} \approx \frac{d\vec{r}}{dt}$ (差分法) |
| **3. 相对速度** | `computeVelBCS` | `src/core/velocity.cpp` | 合成 $\vec{V}_{rel}$ 并投影到叶片坐标系 |
| **4. 尾迹对流** | `AdvanceWakeStructure` | `src/core/wake.cpp` | 欧拉积分更新尾迹点位置 |
| **5. 诱导速度** | `computeInducedVelocity`| `src/core/wake.cpp` | **Biot-Savart 定律**：计算 $\vec{V}_{ind}$ |
| **6. 气动载荷** | `kuttaJoukowskiIteration`| `src/core/wake.cpp` | **Kutta-Joukowski 迭代**：求解 $\Gamma$ |

### 2.2 关键实现细节

#### 2.2.1 运动速度的差分算法
*   **物理对应**：$\vec{V}_{motion} = d\vec{r}/dt$
*   **文件**: `src/core/velocity.cpp` -> `ctdiff`
*   **代码解读**:
    ```cpp
    // 核心逻辑：中心差分 (Central Difference)，利用 t+1 和 t-1 时刻的位置
    if (t >= 1 && t < nTimesteps - 1) {
        result = (pos[t + 1] - pos[t - 1]) * (1.0 / (2.0 * dt));
    }
    // 注意：在 t=N-1 (最后一步) 时，由于没有未来数据，程序自动降级为后向差分。
    // 这可能导致计算出的载荷 (Ct) 在最后一步出现微小跳变，属正常数值现象。
    ```

#### 2.2.2 诱导速度计算
*   **物理对应**：Biot-Savart Law
*   **文件**: `src/core/wake.cpp` -> `computeInducedVelocity`
*   **代码解读**:
    ```cpp
    // 遍历所有涡线段计算诱导速度
    // 在分母中加入 smoothing_term (截断半径) 以防止奇异性 (r -> 0)
    double denominator = r1_r2 * (r1_r2 + dot_r1_r2) + smoothing_term; 
    ```

---

## Part 3: 运行指南 (Run Guide)

### 3.1 编译与运行
```bash
# 1. 编译
mkdir build && cd build
cmake ..
make -j4

# 2. 运行 (默认使用 ../config.json)
./fvw_cpp

# 3. 后处理 (分析结果)
cd ../python
python example.py ../build/results/case_baseline/wake.h5
```

### 3.2 配置文件
主要配置文件位于项目根目录的 `config.json`。常用修改项：
*   `turbine.windSpeed`: 修改风速
*   `simulation.dt`: 修改时间步长
*   `simulation.totalTime`: 修改仿真总时长

更详细的配置说明请参考 [Configuration Guide](configuration.md)。
