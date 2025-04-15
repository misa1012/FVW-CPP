#include "velocity.h"
#include <cmath>

namespace fvw
{

    VelocityData::VelocityData(int nBlades_, int timesteps_, int nShed_)
        : nBlades(nBlades_), nTimesteps(timesteps_), nShed(nShed_)
    {
        bound.resize(nBlades * nTimesteps * nShed);
        blade.resize(nBlades * nTimesteps * nShed);
    }

    const Vec3 &VelocityData::boundAt(int b, int t, int i) const
    {
        return bound[b * nTimesteps * nShed + t * nShed + i];
    }

    const Vec3 &VelocityData::bladeAt(int b, int t, int i) const
    {
        return blade[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &VelocityData::setBoundAt(int b, int t, int i)
    {
        return bound[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &VelocityData::setBladeAt(int b, int t, int i)
    {
        return blade[b * nTimesteps * nShed + t * nShed + i];
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

    // 中心差分
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

    void computeVelocities(VelocityData &vel, NodeAxes &axes, const PositionData &pos,
                           const SimParams &simParams, const TurbineParams &turbineParams)
    {
        // Time array
        std::vector<double> times(simParams.timesteps);
        for (int t = 0; t < simParams.timesteps; ++t)
        {
            times[t] = t * simParams.dt;
        }

        // Compute vel_bound
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
                    vel.setBoundAt(b, t, i) = ctdiff(bound_t, times, t, simParams.timesteps) * -1.0 +
                                              Vec3(turbineParams.windSpeed, 0.0, 0.0);
                }
            }
        }

        // Compute node axes
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

        // Compute vel_blade
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int t = 0; t < simParams.timesteps; ++t)
            {
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    Vec3 &blade = vel.setBladeAt(b, t, i);
                    blade.x = axes.bxnAt(b, t, i).dot(vel.boundAt(b, t, i));
                    blade.y = axes.bynAt(b, t, i).dot(vel.boundAt(b, t, i));
                    blade.z = axes.bznAt(b, t, i).dot(vel.boundAt(b, t, i));
                }
            }
        }
    }

    void computeAoAG(std::vector<std::vector<std::vector<double>>> &aoag, const VelocityData &vel,
                     int nBlades, int timesteps, int nShed)
    {
        for (int b = 0; b < nBlades; ++b)
        {
            for (int t = 0; t < timesteps; ++t)
            {
                for (int i = 0; i < nShed; ++i)
                {
                    aoag[b][t][i] = std::atan2(-vel.bladeAt(b, t, i).y, vel.bladeAt(b, t, i).x) * 180.0 / M_PI;
                    if (std::isnan(aoag[b][t][i]) || std::abs(aoag[b][t][i]) == 180.0)
                    {
                        aoag[b][t][i] = 0.0;
                    }
                }
            }
        }
    }

} // namespace fvw