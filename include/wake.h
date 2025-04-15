#ifndef FVW_WAKE_H
#define FVW_WAKE_H

#include "utils.h"
#include "geometry.h"
#include "position.h"
#include "velocity.h"
#include "bem.h"
#include <vector>

namespace fvw {

struct WakeData {
    // 尾迹涡量：shedding（时间方向），trailing（空间方向）
    std::vector<std::vector<std::vector<double>>> gamma_shed; // [nBlades][nTimesteps][nShed]
    std::vector<std::vector<std::vector<double>>> gamma_trail; // [nBlades][nTimesteps][nShed+1]
    // 尾迹涡点位置
    std::vector<std::vector<std::vector<Vec3>>> wake_pos_shed; // [nBlades][nTimesteps][nShed]
    std::vector<std::vector<std::vector<Vec3>>> wake_pos_trail; // [nBlades][nTimesteps][nShed+1]
    
    WakeData(int nBlades, int nTimesteps, int nShed);
};

// 初始化尾迹（计算初始环量和位置）
void initializeWake(WakeData& wake,
                    const PositionData& pos,
                    const VelocityData& vel,
                    const PerformanceData& perf,
                    const BladeGeometry& geom,
                    const TurbineParams& turbineParams,
                    const SimParams& simParams);

// 更新尾迹位置（对流）
void updateWake(WakeData& wake,
                const VelocityData& vel,
                const SimParams& simParams);

// 计算尾迹诱导速度（Biot-Savart）
void computeInducedVelocity(VelocityData& vel,
                           const WakeData& wake,
                           const BladeGeometry& geom,
                           const TurbineParams& turbineParams,
                           const SimParams& simParams);

} // namespace fvw

#endif // FVW_WAKE_H