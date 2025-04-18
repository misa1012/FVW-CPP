#ifndef FVW_WAKE_H
#define FVW_WAKE_H

#include "geometry.h"
#include "performance.h"
#include "position.h"
#include "utils.h"
#include <vector>

namespace fvw
{

    struct VortexNode
    {
        Vec3 position;
        Vec3 velocity;
    };

    struct VortexLine
    {
        int startNodeIdx;
        int endNodeIdx;
        double gamma;
        bool isShedding; // true: shedding, false: trailing
    };

    struct Wake
    {
        std::vector<std::vector<VortexNode>> nodes; // [timestep][nodes]
        std::vector<std::vector<VortexLine>> lines; // [timestep][lines]
        int nBlades;
        int nShed;
        int nTrail;

        Wake(int nBlades_, int nShed_, int nTrail_)
            : nBlades(nBlades_), nShed(nShed_), nTrail(nTrail_) {}
    };

    void initializeWake(Wake& wake, const BladeGeometry& geom, const PerformanceData& perf,
        const TurbineParams& turbineParams, const PositionData& pos);

    void computeInducedVelocity(std::vector<Vec3> &inducedVel, const Wake &wake,
                                const TurbineParams &turbineParams, int currentTimestep, double cutOff = 0.001);

    // void updateWake(Wake& wake, const PerformanceData& perf, const BladeGeometry& geom,
    //                 const TurbineParams& turbineParams, const SimParams& simParams, int currentTimestep);

} // namespace fvw

#endif // FVW_WAKE_H