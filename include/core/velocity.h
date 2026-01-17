#ifndef FVW_VELOCITY_H
#define FVW_VELOCITY_H

#include "core/utils.h"
#include "core/position.h"
#include <vector>

namespace fvw
{

    class VelICS
    {
    private:
        std::vector<Vec3> data; // ICS velocity at shedding nodes
        int nBlades, nTimesteps, nShed;

    public:
        VelICS(int nBlades_, int timesteps_, int nShed_);
        const Vec3 &at(int b, int t, int i) const;
        Vec3 &setAt(int b, int t, int i);
    };

    class VelBCS
    {
    private:
        std::vector<Vec3> data; // BCS velocity at shedding nodes
        int nBlades, nTimesteps, nShed;

    public:
        VelBCS(int nBlades_, int timesteps_, int nShed_);
        const Vec3 &at(int b, int t, int i) const;
        Vec3 &setAt(int b, int t, int i);
    };

    class NodeAxes
    {
    private:
        std::vector<Vec3> bxn, byn, bzn; // Shedding node axes
        std::vector<Vec3> bxt, byt, bzt; // Trailing node axes
        int nBlades, nTimesteps, nTrail, nShed;

    public:
        NodeAxes(int nBlades_, int timesteps_, int nTrail_, int nShed_);
        Vec3 &bxnAt(int b, int t, int i);
        Vec3 &bynAt(int b, int t, int i);
        Vec3 &bznAt(int b, int t, int i);
        Vec3 &bxtAt(int b, int t, int i);
        Vec3 &bytAt(int b, int t, int i);
        Vec3 &bztAt(int b, int t, int i);
    };

    void computeVelICS(VelICS &velICS, const PositionData &pos,
                       const SimParams &simParams, const TurbineParams &turbineParams);

    void computeVelBCS(VelBCS &velBCS, const VelICS &velICS, NodeAxes &axes, const PositionData &pos,
                       const SimParams &simParams, const TurbineParams &turbineParams);

} // namespace fvw

#endif // FVW_VELOCITY_H