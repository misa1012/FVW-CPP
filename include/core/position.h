#ifndef FVW_POSITION_H
#define FVW_POSITION_H

#include "core/utils.h"
#include "core/geometry.h"
#include <vector>
#include <string>

namespace fvw
{

    class PositionData
    {
    private:
        std::vector<Vec3> lead, quarter, trail, colloc, bound, end;
        std::vector<Vec3> hub, platform;
        int nBlades, nTimesteps, nTrail, nShed;

    public:
        PositionData(int nBlades_, int timesteps_, int nTrail_, int nShed_);

        // 只读访问器
        const Vec3 &leadAt(int b, int t, int i) const;
        const Vec3 &quarterAt(int b, int t, int i) const;
        const Vec3 &trailAt(int b, int t, int i) const;
        const Vec3 &collocAt(int b, int t, int i) const;
        const Vec3 &boundAt(int b, int t, int i) const;
        const Vec3 &endAt(int b, int t, int i) const;
        const Vec3 &hubAt(int t) const;
        const Vec3 &platformAt(int t) const;

        // 可写访问器
        Vec3 &setLeadAt(int b, int t, int i);
        Vec3 &setQuarterAt(int b, int t, int i);
        Vec3 &setTrailAt(int b, int t, int i);
        Vec3 &setCollocAt(int b, int t, int i);
        Vec3 &setBoundAt(int b, int t, int i);
        Vec3 &setEndAt(int b, int t, int i);
        Vec3 &setHubAt(int t);
        Vec3 &setPlatformAt(int t);
    };

    // 方向余弦矩阵旋转
    Vec3 DCMRot(const Vec3 &x, const std::vector<double> &t,
                const std::vector<double> &A_init, const std::string &rotseq, int rev);

    // 计算位置
    void computePositions(PositionData &pos, const SimParams &simParams,
                          const TurbineParams &turbineParams, const BladeGeometry &geom);

} // namespace fvw

#endif // FVW_POSITION_H