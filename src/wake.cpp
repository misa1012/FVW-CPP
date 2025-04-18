#include "wake.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

namespace fvw
{
    // The first wake
    void initializeWake(Wake &wake, const BladeGeometry &geom, const PerformanceData &perf,
                        const TurbineParams &turbineParams, const PositionData &pos)
    {
        wake.nodes.clear();
        wake.lines.clear();
        wake.nodes.emplace_back();
        wake.lines.emplace_back();

        std::vector<VortexNode> &nodes = wake.nodes[0];
        std::vector<VortexLine> &lines = wake.lines[0];

        nodes.reserve(wake.nBlades * 2 * wake.nTrail);
        lines.reserve(wake.nBlades * (2 * wake.nShed + wake.nTrail));

        std::vector<std::vector<int>> controlNodeIdx(wake.nBlades, std::vector<int>(wake.nTrail));
        std::vector<std::vector<int>> trailNodeIdx(wake.nBlades, std::vector<int>(wake.nTrail));

        // 添加节点（使用全局坐标）
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nTrail; ++i)
            {
                controlNodeIdx[b][i] = nodes.size();
                nodes.push_back({pos.quarterAt(b, 0, i), Vec3(turbineParams.windSpeed, 0.0, 0.0)});
                trailNodeIdx[b][i] = nodes.size();
                nodes.push_back({pos.trailAt(b, 0, i), Vec3(turbineParams.windSpeed, 0.0, 0.0)});
            }
        }

        // Bound vortex: n_shed lines per blade,首尾连接
        std::vector<std::vector<double>> gamma_shed(wake.nBlades, std::vector<double>(wake.nShed));
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nShed; ++i)
            {
                double cl = perf.clAt(b, 0, i);
                double chord = geom.chordShedding[i];
                double V_rel = turbineParams.windSpeed;
                gamma_shed[b][i] = 0.5 * V_rel * chord * cl;

                int startIdx = controlNodeIdx[b][i];
                int endIdx = controlNodeIdx[b][i + 1]; // n_trail = n_shed + 1
                lines.push_back({startIdx, endIdx, gamma_shed[b][i], true});
            }
        }

        // Trailing vortex: n_trail lines per blade
        std::vector<std::vector<double>> gamma_trail(wake.nBlades, std::vector<double>(wake.nTrail));
        for (int b = 0; b < wake.nBlades; ++b)
        {
            // 边界填充 0
            gamma_trail[b][0] = 0.0;
            for (int i = 0; i < wake.nShed; ++i)
            {
                gamma_trail[b][i + 1] = gamma_shed[b][i];
            }
            for (int i = 1; i < wake.nTrail; ++i)
            {
                gamma_trail[b][i] = gamma_trail[b][i] - gamma_trail[b][i - 1];
            }
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = controlNodeIdx[b][i];
                int endIdx = trailNodeIdx[b][i];
                lines.push_back({startIdx, endIdx, gamma_trail[b][i], false});
            }
        }

        // Shedding vortex: n_shed lines per blade
        std::vector<std::vector<double>> gamma(wake.nBlades, std::vector<double>(wake.nShed));
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = trailNodeIdx[b][i];
                int endIdx = trailNodeIdx[b][i + 1];
                lines.push_back({startIdx, endIdx, -1 * gamma_shed[b][i], true});
            }
        }

        // 计算诱导速度，叠加到节点速度
        std::vector<Vec3> inducedVel(nodes.size());
        computeInducedVelocity(inducedVel, wake, turbineParams, 0);
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            nodes[i].velocity = nodes[i].velocity + inducedVel[i];
        }

        std::cout << "Wake initialized: nBlades=" << wake.nBlades
                  << ", nShed=" << wake.nShed
                  << ", nTrail=" << wake.nTrail
                  << ", nodes[0]=" << nodes.size()
                  << ", lines[0]=" << lines.size() << std::endl
                  << "First node vel=" << to_string(nodes[0].velocity) << std::endl;
    }

    void computeInducedVelocity(std::vector<Vec3> &inducedVel, const Wake &wake,
                                const TurbineParams &turbineParams, int currentTimestep, double cutOff)
    {
        const auto &nodes = wake.nodes[currentTimestep];
        const auto &lines = wake.lines[currentTimestep];
        inducedVel.resize(nodes.size(), Vec3(0.0, 0.0, 0.0));

        double cutOffScaled = cutOff * turbineParams.rTip;
        double denominatorThreshold = 1e-10; // 奇异点阈值

        for (size_t n = 0; n < nodes.size(); ++n)
        {
            Vec3 p = nodes[n].position;
            Vec3 vel(0.0, 0.0, 0.0);

            for (size_t l = 0; l < lines.size(); ++l)
            {
                const VortexLine &line = lines[l];
                Vec3 x1 = nodes[line.startNodeIdx].position;
                Vec3 x2 = nodes[line.endNodeIdx].position;
                double gamma = line.gamma;

                Vec3 l_vec = x2 - x1; // 改为 l_vec，避免与循环变量 l 冲突
                double l_squared = l_vec.x * l_vec.x + l_vec.y * l_vec.y + l_vec.z * l_vec.z;
                double cut_l = cutOffScaled * cutOffScaled * l_squared;
                double coeff = gamma / (4.0 * M_PI);

                Vec3 r1 = p - x1;
                Vec3 r2 = p - x2;
                double r1_norm = r1.norm();
                double r2_norm = r2.norm();
                double r1_r2 = r1_norm * r2_norm;
                double dot_r1_r2 = r1.dot(r2);
                Vec3 cross_r1_r2 = r1.cross(r2);

                double denominator = r1_r2 * (r1_r2 + dot_r1_r2) + cut_l;

                // 跳过奇异点
                if (denominator < denominatorThreshold || r1_norm < 1e-6 || r2_norm < 1e-6)
                {
                    if (n == 0 && l < 5)
                    {
                        std::cout << "Line " << l << ": Skipped due to small denominator=" << denominator
                                  << ", r1_norm=" << r1_norm << ", r2_norm=" << r2_norm << std::endl;
                    }
                    continue;
                }

                double contribution = coeff * (r1_norm + r2_norm) / denominator;
                Vec3 vel_contrib = cross_r1_r2 * contribution;
                vel = vel + cross_r1_r2 * contribution;
            }

            inducedVel[n] = vel;

            std::cout << "Induced velocity at t=" << currentTimestep
                      << ": vel[" << n << "]=" << to_string(inducedVel[n]) << std::endl;
        }
    }
} // namespace fvw