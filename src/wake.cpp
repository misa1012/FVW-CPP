#include "wake.h"
#include <cmath>

namespace fvw
{

    WakeData::WakeData(int nBlades, int nTimesteps, int nShed)
    {
        gamma_shed.resize(nBlades, std::vector<std::vector<double>>(
                                       nTimesteps, std::vector<double>(nShed, 0.0)));
        gamma_trail.resize(nBlades, std::vector<std::vector<double>>(
                                        nTimesteps, std::vector<double>(nShed + 1, 0.0)));
        wake_pos_shed.resize(nBlades, std::vector<std::vector<Vec3>>(
                                          nTimesteps, std::vector<Vec3>(nShed, Vec3(0.0, 0.0, 0.0))));
        wake_pos_trail.resize(nBlades, std::vector<std::vector<Vec3>>(
                                           nTimesteps, std::vector<Vec3>(nShed + 1, Vec3(0.0, 0.0, 0.0))));
    }

    void initializeWake(WakeData &wake,
                        const PositionData &pos,
                        const VelocityData &vel,
                        const PerformanceData &perf,
                        const BladeGeometry &geom,
                        const TurbineParams &turbineParams,
                        const SimParams &simParams)
    {
        const double pi = M_PI;
        const double rho = turbineParams.rho;

        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int t = 0; t < simParams.timesteps; ++t)
            {
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    double cl = perf.clAt(b, t, i);
                    Vec3 V_rel = vel.bladeAt(b, t, i);
                    double V_mag = sqrt(V_rel.x * V_rel.x + V_rel.y * V_rel.y);
                    double chord = geom.chordShedding[i];
                    wake.gamma_shed[b][t][i] = 0.5 * cl * V_mag * chord / rho;

                    wake.wake_pos_shed[b][t][i] = pos.boundAt(b, t, i);
                }
                for (int i = 0; i <= turbineParams.nSegments; ++i)
                {
                    if (i == 0)
                    {
                        wake.gamma_trail[b][t][i] = -wake.gamma_shed[b][t][i];
                    }
                    else if (i == turbineParams.nSegments)
                    {
                        wake.gamma_trail[b][t][i] = wake.gamma_shed[b][t][i - 1];
                    }
                    else
                    {
                        wake.gamma_trail[b][t][i] = wake.gamma_shed[b][t][i - 1] - wake.gamma_shed[b][t][i];
                    }
                    wake.wake_pos_trail[b][t][i] = pos.trailAt(b, t, i);
                }
            }
        }
    }

    void updateWake(WakeData &wake,
                    const VelocityData &vel,
                    const SimParams &simParams)
    {
        for (int b = 0; b < wake.gamma_shed.size(); ++b)
        {
            for (int t = 1; t < wake.gamma_shed[0].size(); ++t)
            {
                for (int i = 0; i < wake.gamma_shed[0][0].size(); ++i)
                {
                    Vec3 v = vel.boundAt(b, t - 1, i);
                    wake.wake_pos_shed[b][t][i] = wake.wake_pos_shed[b][t - 1][i] + v * simParams.dt;
                }
                for (int i = 0; i < wake.gamma_trail[0][0].size(); ++i)
                {
                    int shed_idx = (i < wake.gamma_shed[0][0].size()) ? i : wake.gamma_shed[0][0].size() - 1;
                    Vec3 v = vel.boundAt(b, t - 1, shed_idx);
                    wake.wake_pos_trail[b][t][i] = wake.wake_pos_trail[b][t - 1][i] + v * simParams.dt;
                }
            }
        }
    }

    void computeInducedVelocity(VelocityData &vel,
                                const WakeData &wake,
                                const BladeGeometry &geom,
                                const TurbineParams &turbineParams,
                                const SimParams &simParams)
    {
        const double pi = M_PI;
        const double epsilon = 1e-6;

        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int t = 0; t < simParams.timesteps; ++t)
            {
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    Vec3 induced_vel(0.0, 0.0, 0.0);
                    Vec3 x = geom.colloc[i];

                    // shedding 涡量
                    for (int wb = 0; wb < wake.gamma_shed.size(); ++wb)
                    {
                        for (int wt = 0; wt < wake.gamma_shed[0].size() - 1; ++wt)
                        {
                            for (int wi = 0; wi < wake.gamma_shed[0][0].size(); ++wi)
                            {
                                Vec3 x1 = wake.wake_pos_shed[wb][wt][wi];
                                Vec3 x2 = wake.wake_pos_shed[wb][wt + 1][wi];
                                double gamma = wake.gamma_shed[wb][wt][wi];
                                Vec3 dl = x2 - x1;
                                Vec3 r = x - x1;
                                double r_norm = r.norm();
                                if (r_norm > epsilon)
                                {
                                    double coef = gamma / (4 * pi * r_norm * r_norm * r_norm);
                                    induced_vel = induced_vel + (r.cross(dl) * coef); // 使用 operator+
                                }
                            }
                        }
                    }

                    // trailing 涡量
                    for (int wb = 0; wb < wake.gamma_trail.size(); ++wb)
                    {
                        for (int wt = 0; wt < wake.gamma_trail[0].size(); ++wt)
                        {
                            for (int wi = 0; wi < wake.gamma_trail[0][0].size() - 1; ++wi)
                            {
                                Vec3 x1 = wake.wake_pos_trail[wb][wt][wi];
                                Vec3 x2 = wake.wake_pos_trail[wb][wt][wi + 1];
                                double gamma = wake.gamma_trail[wb][wt][wi];
                                Vec3 dl = x2 - x1;
                                Vec3 r = x - x1;
                                double r_norm = r.norm();
                                if (r_norm > epsilon)
                                {
                                    double coef = gamma / (4 * pi * r_norm * r_norm * r_norm);
                                    induced_vel = induced_vel + (r.cross(dl) * coef); // 使用 operator+
                                }
                            }
                        }
                    }

                    Vec3 current_vel = vel.boundAt(b, t, i);
                    vel.setBoundAt(b, t, i) = current_vel + induced_vel;
                }
            }
        }
    }

} // namespace fvw