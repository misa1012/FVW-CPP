#ifndef FVW_WAKE_H
#define FVW_WAKE_H

#include "geometry.h"
#include "performance.h"
#include "position.h"
#include "utils.h"
#include "velocity.h"
#include "airfoil.h"
#include <vector>

namespace fvw
{

    struct VortexNode
    {
        Vec3 position;
        Vec3 velocity;
        int idx;
    };

    enum class VortexLineType
    {
        Bound,
        Shed,
        Trailing,
        NewShed
    };

    struct VortexLine
    {
        int startNodeIdx;
        int endNodeIdx;
        double gamma;
        VortexLineType type;
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

    void computeInducedVelocity(std::vector<Vec3> &inducedVel, std::vector<VortexNode> &pNodes, const std::vector<VortexNode> &nodes,
                                const std::vector<VortexLine> &lines, const TurbineParams &turbineParams, double cutOff = 0.001);

    void initializeWake(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                        const TurbineParams &turbineParams, const PositionData &pos, double dt);

    void kuttaJoukowskiIteration(Wake &wake, PerformanceData &perf, const BladeGeometry &geom, NodeAxes &axes,
                                 const TurbineParams &turbineParams, const PositionData &pos, VelBCS &velBCS, std::vector<AirfoilData> &airfoils);

} // namespace fvw

#endif // FVW_WAKE_H