#include "velocity.h"
#include <cmath>

namespace fvw
{

    VelICS::VelICS(int nBlades_, int timesteps_, int nShed_)
        : nBlades(nBlades_), nTimesteps(timesteps_), nShed(nShed_)
    {
        data.resize(nBlades * nTimesteps * nShed);
    }

    const Vec3 &VelICS::at(int b, int t, int i) const
    {
        return data[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &VelICS::setAt(int b, int t, int i)
    {
        return data[b * nTimesteps * nShed + t * nShed + i];
    }

    VelBCS::VelBCS(int nBlades_, int timesteps_, int nShed_)
        : nBlades(nBlades_), nTimesteps(timesteps_), nShed(nShed_)
    {
        data.resize(nBlades * nTimesteps * nShed);
    }

    const Vec3 &VelBCS::at(int b, int t, int i) const
    {
        return data[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &VelBCS::setAt(int b, int t, int i)
    {
        return data[b * nTimesteps * nShed + t * nShed + i];
    }

    NodeAxes::NodeAxes(int nBlades_, int timesteps_, int nTrail_, int nShed_)
        : nBlades(nBlades_), nTimesteps(timesteps_), nTrail(nTrail_), nShed(nShed_)
    {
        bxn.resize(nBlades * nTimesteps * nShed);
        byn.resize(nBlades * nTimesteps * nShed);
        bzn.resize(nBlades * nTimesteps * nShed);
        bxt.resize(nBlades * nTimesteps * nTrail);
        byt.resize(nBlades * nTimesteps * nTrail);
        bzt.resize(nBlades * nTimesteps * nTrail);
    }

    Vec3 &NodeAxes::bxnAt(int b, int t, int i) { return bxn[b * nTimesteps * nShed + t * nShed + i]; }
    Vec3 &NodeAxes::bynAt(int b, int t, int i) { return byn[b * nTimesteps * nShed + t * nShed + i]; }
    Vec3 &NodeAxes::bznAt(int b, int t, int i) { return bzn[b * nTimesteps * nShed + t * nShed + i]; }
    Vec3 &NodeAxes::bxtAt(int b, int t, int i) { return bxt[b * nTimesteps * nTrail + t * nTrail + i]; }
    Vec3 &NodeAxes::bytAt(int b, int t, int i) { return byt[b * nTimesteps * nTrail + t * nTrail + i]; }
    Vec3 &NodeAxes::bztAt(int b, int t, int i) { return bzt[b * nTimesteps * nTrail + t * nTrail + i]; }

    static Vec3 ctdiff(const std::vector<Vec3> &pos, const std::vector<double> &times, int t, int nTimesteps)
    {
        if (t == 0)
        {
            return (pos[1] - pos[0]) * (1.0 / (times[1] - times[0]));
        }
        else if (t == nTimesteps - 1)
        {
            return (pos[t] - pos[t - 1]) * (1.0 / (times[t] - times[t - 1]));
        }
        else
        {
            return (pos[t + 1] - pos[t - 1]) * (0.5 / (times[t + 1] - times[t - 1]));
        }
    }

    void computeVelICS(VelICS &velICS, const PositionData &pos,
                       const SimParams &simParams, const TurbineParams &turbineParams)
    {
        // Compute time array
        std::vector<double> times(simParams.timesteps);
        for (int t = 0; t < simParams.timesteps; ++t)
        {
            times[t] = t * simParams.dt;
        }

        // Compute ICS velocities using center difference and wind speed
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int i = 0; i < turbineParams.nSegments; ++i)
            {
                std::vector<Vec3> bound_t(simParams.timesteps);
                for (int t = 0; t < simParams.timesteps; ++t)
                {
                    bound_t[t] = pos.boundAt(b, t, i);
                }
                for (int t = 0; t < simParams.timesteps; ++t)
                {
                    velICS.setAt(b, t, i) = ctdiff(bound_t, times, t, simParams.timesteps) * -1.0 +
                                            Vec3(turbineParams.windSpeed, 0.0, 0.0);
                }
            }
        }
    }

    void computeVelBCS(VelBCS &velBCS, const VelICS &velICS, NodeAxes &axes, const PositionData &pos,
                       const SimParams &simParams, const TurbineParams &turbineParams)
    {
        // Compute node axes for trailing and shedding nodes
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int t = 0; t < simParams.timesteps; ++t)
            {
                // Trailing nodes
                for (int i = 0; i < turbineParams.nSegments + 1; ++i)
                {
                    axes.bxtAt(b, t, i) = pos.trailAt(b, t, i) - pos.quarterAt(b, t, i);
                    double norm = axes.bxtAt(b, t, i).norm();
                    if (norm > 1e-10)
                        axes.bxtAt(b, t, i) = axes.bxtAt(b, t, i) * (1.0 / norm);
                }
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    axes.bztAt(b, t, i) = pos.quarterAt(b, t, i + 1) - pos.quarterAt(b, t, i);
                    double norm = axes.bztAt(b, t, i).norm();
                    if (norm > 1e-10)
                        axes.bztAt(b, t, i) = axes.bztAt(b, t, i) * (1.0 / norm);
                }
                axes.bztAt(b, t, turbineParams.nSegments) = axes.bztAt(b, t, turbineParams.nSegments - 1);
                for (int i = 0; i < turbineParams.nSegments + 1; ++i)
                {
                    axes.bytAt(b, t, i) = axes.bztAt(b, t, i).cross(axes.bxtAt(b, t, i));
                    double norm = axes.bytAt(b, t, i).norm();
                    if (norm > 1e-10)
                        axes.bytAt(b, t, i) = axes.bytAt(b, t, i) * (1.0 / norm);
                }

                // Shedding nodes
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    axes.bxnAt(b, t, i) = pos.endAt(b, t, i) - pos.boundAt(b, t, i);
                    double norm = axes.bxnAt(b, t, i).norm();
                    if (norm > 1e-10)
                        axes.bxnAt(b, t, i) = axes.bxnAt(b, t, i) * (1.0 / norm);
                }
                for (int i = 0; i < turbineParams.nSegments - 1; ++i)
                {
                    axes.bznAt(b, t, i) = pos.boundAt(b, t, i + 1) - pos.boundAt(b, t, i);
                    double norm = axes.bznAt(b, t, i).norm();
                    if (norm > 1e-10)
                        axes.bznAt(b, t, i) = axes.bznAt(b, t, i) * (1.0 / norm);
                }
                axes.bznAt(b, t, turbineParams.nSegments - 1) = axes.bznAt(b, t, turbineParams.nSegments - 2);
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    axes.bynAt(b, t, i) = axes.bznAt(b, t, i).cross(axes.bxnAt(b, t, i));
                    double norm = axes.bynAt(b, t, i).norm();
                    if (norm > 1e-10)
                        axes.bynAt(b, t, i) = axes.bynAt(b, t, i) * (1.0 / norm);
                }
            }
        }

        // Compute BCS velocities by projecting ICS velocities onto node axes
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int t = 0; t < simParams.timesteps; ++t)
            {
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    Vec3 &bcs = velBCS.setAt(b, t, i);
                    bcs.x = axes.bxnAt(b, t, i).dot(velICS.at(b, t, i));
                    bcs.y = axes.bynAt(b, t, i).dot(velICS.at(b, t, i));
                    bcs.z = axes.bznAt(b, t, i).dot(velICS.at(b, t, i));
                }
            }
        }
    }

} // namespace fvw