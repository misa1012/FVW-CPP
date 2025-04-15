#ifndef FVW_VELOCITY_H
#define FVW_VELOCITY_H

#include "utils.h"
#include "position.h"
#include <vector>

namespace fvw {

class VelocityData {
private:
    std::vector<Vec3> bound; // ICS velocity at shedding nodes
    std::vector<Vec3> blade; // BCS velocity at shedding nodes
    int nBlades, nTimesteps, nShed;

public:
    VelocityData(int nBlades_, int timesteps_, int nShed_);

    // 只读访问器
    const Vec3& boundAt(int b, int t, int i) const;
    const Vec3& bladeAt(int b, int t, int i) const;

    // 可写访问器
    Vec3& setBoundAt(int b, int t, int i);
    Vec3& setBladeAt(int b, int t, int i);
};

class NodeAxes {
private:
    std::vector<Vec3> bxn, byn, bzn; // Shedding node axes
    std::vector<Vec3> bxt, byt, bzt; // Trailing node axes
    int nBlades, nTimesteps, nTrail, nShed;

public:
    NodeAxes(int nBlades_, int timesteps_, int nTrail_, int nShed_);
    Vec3& bxnAt(int b, int t, int i);
    Vec3& bynAt(int b, int t, int i);
    Vec3& bznAt(int b, int t, int i);
    Vec3& bxtAt(int b, int t, int i);
    Vec3& bytAt(int b, int t, int i);
    Vec3& bztAt(int b, int t, int i);
};

// 计算速度和节点基向量
void computeVelocities(VelocityData& vel, NodeAxes& axes, const PositionData& pos,
                       const SimParams& simParams, const TurbineParams& turbineParams);

// 计算几何迎角
void computeAoAG(std::vector<std::vector<std::vector<double>>>& aoag, const VelocityData& vel,
                 int nBlades, int timesteps, int nShed);

} // namespace fvw

#endif // FVW_VELOCITY_H