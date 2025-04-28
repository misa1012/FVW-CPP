#include "wake.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

namespace fvw
{
    // Biot-Savart function
    void computeInducedVelocity(std::vector<Vec3> &inducedVel, const std::vector<VortexNode> &nodes,
                                const std::vector<VortexLine> &lines, const TurbineParams &turbineParams, double cutOff)
    {
        inducedVel.resize(nodes.size(), Vec3(0.0, 0.0, 0.0));

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

                Vec3 l_vec = x2 - x1;
                double l_squared = l_vec.x * l_vec.x + l_vec.y * l_vec.y + l_vec.z * l_vec.z;
                double cut_l = cutOff * cutOff * l_squared;
                double coeff = gamma / (4.0 * M_PI);

                Vec3 r1 = p - x1;
                Vec3 r2 = p - x2;
                double r1_norm = r1.norm();
                double r2_norm = r2.norm();
                double r1_r2 = r1_norm * r2_norm;
                double dot_r1_r2 = r1.dot(r2);
                Vec3 cross_r1_r2 = r1.cross(r2);

                double denominator = r1_r2 * (r1_r2 + dot_r1_r2) + cut_l;

                double contribution = coeff * (r1_norm + r2_norm) / denominator;
                Vec3 vel_contrib = cross_r1_r2 * contribution;

                vel = vel + vel_contrib;
            }

            inducedVel[n] = vel;

            // std::cout << "Induced velocity at t=" << currentTimestep
            //           << ": vel[" << n << "]=" << to_string(inducedVel[n]) << std::endl;
        }
    }

    // The first wake
    void initializeWake(Wake &wake, const BladeGeometry &geom, const PerformanceData &perf,
                        const TurbineParams &turbineParams, const PositionData &pos, double dt)
    {
        bool if_verbose = true;

        wake.nodes.clear();
        wake.lines.clear();
        wake.nodes.emplace_back();
        wake.lines.emplace_back();

        std::vector<VortexNode> &nodes = wake.nodes[0];
        std::vector<VortexLine> &lines = wake.lines[0];

        int NumNode0 = wake.nBlades * 2 * wake.nTrail;
        int NumLine0 = wake.nBlades * (2 * wake.nShed + wake.nTrail);

        nodes.reserve(NumNode0);
        lines.reserve(NumLine0);

        std::vector<std::vector<int>> boundNodeIdx(wake.nBlades, std::vector<int>(wake.nTrail));
        std::vector<std::vector<int>> trailNodeIdx(wake.nBlades, std::vector<int>(wake.nTrail));

        // Initialize wake node
        int node_idx_counter = 0;
        Vec3 initialVel = Vec3(turbineParams.windSpeed, 0.0, 0.0);
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nTrail; ++i)
            {
                boundNodeIdx[b][i] = node_idx_counter;
                nodes.push_back({pos.quarterAt(b, 0, i), initialVel, node_idx_counter++});
                trailNodeIdx[b][i] = node_idx_counter;
                nodes.push_back({pos.trailAt(b, 0, i), initialVel, node_idx_counter++});

                // if (if_verbose && b == 0) {
                //     std::cout << "i=" << std::setw(2) << i
                //               << ", controlNodeIdx=" << std::setw(3) << controlNodeIdx[b][i]
                //               << ", pos=" << to_string(nodes[controlNodeIdx[b][i]].position)
                //               << ", vel=" << to_string(nodes[controlNodeIdx[b][i]].velocity) << std::endl;
                //     std::cout << "i=" << std::setw(2) << i
                //               << ", trailNodeIdx=" << std::setw(3) << trailNodeIdx[b][i]
                //               << ", pos=" << to_string(nodes[trailNodeIdx[b][i]].position)
                //               << ", vel=" << to_string(nodes[trailNodeIdx[b][i]].velocity) << std::endl;
                // }
            }
        }

        // Initialize wake line
        // Bound vortex: n_shed lines per blade,首尾连接
        std::vector<std::vector<double>> gamma_bound(wake.nBlades, std::vector<double>(wake.nShed));
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nShed; ++i)
            {
                double cl = perf.clAt(b, 0, i);
                double chord = geom.chordShedding[i];
                gamma_bound[b][i] = 0.5 * turbineParams.windSpeed * chord * cl;

                int startIdx = boundNodeIdx[b][i];
                int endIdx = boundNodeIdx[b][i + 1]; // n_trail = n_shed + 1
                lines.push_back({startIdx, endIdx, gamma_bound[b][i], VortexLineType::Bound});
            }
        }

        // Trailing vortex: n_trail lines per blade
        std::vector<std::vector<double>> gamma_trail(wake.nBlades, std::vector<double>(wake.nTrail));
        for (int b = 0; b < wake.nBlades; ++b)
        {
            gamma_trail[b][0] = gamma_bound[b][0];
            for (int i = 1; i < wake.nShed; ++i)
            {
                gamma_trail[b][i] = gamma_bound[b][i] - gamma_bound[b][i - 1];
            }
            gamma_trail[b][wake.nShed] = -gamma_bound[b][wake.nShed - 1];

            // 这里要注意，debug过，要先trail再control，这会影响trailing induced velocity的正负号
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = trailNodeIdx[b][i];
                int endIdx = boundNodeIdx[b][i];
                lines.push_back({startIdx, endIdx, gamma_trail[b][i], VortexLineType::Trailing});
            }
        }

        // Shedding vortex: n_shed lines per blade
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = trailNodeIdx[b][i];
                int endIdx = trailNodeIdx[b][i + 1];
                lines.push_back({startIdx, endIdx, -1 * gamma_bound[b][i], VortexLineType::Shed});
            }
        }

        if (if_verbose)
        {
            std::cout << "Wake initialized: nBlades=" << wake.nBlades
                      << ", nShed=" << wake.nShed
                      << ", nTrail=" << wake.nTrail
                      << ", nodes[0]=" << nodes.size()
                      << ", lines[0]=" << lines.size() << std::endl;
        }

        // if(if_verbose)// 独立的涡量验证输出
        // {
        //     std::cout << std::fixed << std::setprecision(6);
        //     // 输出 gamma_bound
        //     std::cout << "[Gamma] Blade 0, t=0, gamma_bound:" << std::endl;
        //     for (int i = 0; i < wake.nShed; ++i) {
        //         double cl = perf.clAt(0, 0, i);
        //         double chord = geom.chordShedding[i];
        //         double V_rel = turbineParams.windSpeed;
        //         double gamma = gamma_bound[0][i];
        //         if (std::isnan(gamma) || std::abs(cl) > 10.0) {
        //             std::cerr << "[Warning] Invalid gamma_bound at i=" << i
        //                       << ", cl=" << cl << ", gamma=" << gamma << std::endl;
        //         }
        //         std::cout << "i=" << std::setw(2) << i
        //                   << ", cl=" << std::setw(10) << cl
        //                   << ", chord=" << std::setw(10) << chord
        //                   << ", V_rel=" << std::setw(10) << V_rel
        //                   << ", gamma_bound=" << std::setw(10) << gamma << std::endl;
        //     }

        //     // 输出 gamma_trail
        //     std::cout << "[Gamma] Blade 0, t=0, gamma_trail:" << std::endl;
        //     for (int i = 0; i < wake.nTrail; ++i) {
        //         double gamma = gamma_trail[0][i];
        //         if (std::isnan(gamma)) {
        //             std::cerr << "[Warning] Invalid gamma_trail at i=" << i
        //                       << ", gamma=" << gamma << std::endl;
        //         }
        //         std::cout << "i=" << std::setw(2) << i
        //                   << ", gamma_trail=" << std::setw(10) << gamma << std::endl;
        //     }

        //     // 输出 gamma_shedding
        //     std::cout << "[Gamma] Blade 0, t=0, gamma_shedding:" << std::endl;
        //     for (int i = 0; i < wake.nShed; ++i) {
        //         double gamma = -1 * gamma_bound[0][i];
        //         if (std::isnan(gamma)) {
        //             std::cerr << "[Warning] Invalid gamma_shedding at i=" << i
        //                       << ", gamma=" << gamma << std::endl;
        //         }
        //         std::cout << "i=" << std::setw(2) << i
        //                   << ", gamma_shedding=" << std::setw(10) << gamma << std::endl;
        //     }
        // }

        // 计算诱导速度，叠加到节点速度
        // Update velocity
        std::vector<Vec3> inducedVel(nodes.size());
        std::vector<Vec3> newTrailPos(nodes.size());

        computeInducedVelocity(inducedVel, wake.nodes[0], wake.lines[0], turbineParams);

        for (int i = 0; i < NumNode0; ++i)
        {
            nodes[i].velocity = nodes[i].velocity + inducedVel[i];
            newTrailPos[i] = nodes[i].position + nodes[i].velocity * dt;
        }

        // if (if_verbose)
        // {
        //     int blade = 0;
        //     int timestep = 0;
        //     // 输出指定叶片和时间步的诱导速度 x 坐标
        //     std::cout << std::fixed << std::setprecision(6);
        //     std::cout << "[Induced Velocity] Blade " << blade << ", t=" << timestep
        //               << ", inducedVel.x:" << std::endl;
        //     for (int i = 0; i < wake.nTrail; ++i)
        //     {
        //         int controlIdx = controlNodeIdx[blade][i];
        //         int trailIdx = trailNodeIdx[blade][i];
        //         std::cout << "i=" << std::setw(2) << i
        //                   << ", controlNodeIdx=" << std::setw(3) << controlIdx
        //                   << ", inducedVel.x=" << std::setw(10) << inducedVel[controlIdx].x << std::endl;
        //         std::cout << "i=" << std::setw(2) << i
        //                   << ", trailNodeIdx=" << std::setw(3) << trailIdx
        //                   << ", inducedVel.x=" << std::setw(10) << inducedVel[trailIdx].x << std::endl;
        //     }

        //     // 输出指定叶片的第一个节点的 velocity
        //     int firstNodeIdx = controlNodeIdx[blade][0];
        //     std::cout << "First node vel (Blade " << blade << ", t=" << timestep
        //               << ")=" << to_string(nodes[firstNodeIdx].velocity) << std::endl;
        // }

        // 第一步对流：使用前向欧拉法更新 trailing vortex 节点
        std::cout << "Advancing timestep = 1" << std::endl;

        wake.nodes.emplace_back(); // 创建 t=1 的节点数组，这之后不能再引用lines了，因为很可能对wake已经重新分配了内存
        std::vector<VortexNode> &nodesNext = wake.nodes[1];
        int NumNode1 = wake.nBlades * 3 * wake.nTrail;
        nodesNext.reserve(NumNode1);

        for (int i = 0; i < NumNode0; ++i)
        {
            nodesNext.push_back({newTrailPos[i], initialVel, wake.nodes[0][i].idx});
        }

        // 对于之前的Node：
        // Convect all nodes forward in time using Forward Euler
        // 对流所有 t=0 节点（绑定和尾迹节点）并添加新绑定节点

        std::vector<std::vector<int>> newBoundNodeIdx(wake.nBlades, std::vector<int>(wake.nTrail));
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nTrail; ++i)
            {
                Vec3 newBoundPosT1 = pos.quarterAt(b, 1, i); // t=1 的 1/4 弦长位置
                newBoundNodeIdx[b][i] = node_idx_counter;
                nodesNext.push_back({newBoundPosT1, initialVel, node_idx_counter++});
            }
        }

        // Update line
        wake.lines.emplace_back();
        std::vector<VortexLine> &linesNext = wake.lines[1];
        int NumLine1 = wake.nBlades * (3 * wake.nShed + 2 * wake.nTrail);
        linesNext.reserve(NumLine1);

        for (int n = 0; n < NumLine0; ++n)
        {
            if (wake.lines[0][n].type != VortexLineType::Bound)
            {
                linesNext.push_back(wake.lines[0][n]);
            }
        }

        // Update bound vortex
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = newBoundNodeIdx[b][i];
                int endIdx = newBoundNodeIdx[b][i + 1];
                linesNext.push_back({startIdx, endIdx, gamma_bound[b][i], VortexLineType::Bound});
            }
        }

        // Update trail vortex
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = boundNodeIdx[b][i];
                int endIdx = newBoundNodeIdx[b][i];
                linesNext.push_back({startIdx, endIdx, gamma_trail[b][i], VortexLineType::Trailing});
            }
        }

        // Update shed vortex
        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = boundNodeIdx[b][i];
                int endIdx = boundNodeIdx[b][i + 1];
                linesNext.push_back({startIdx, endIdx, -1 * gamma_bound[b][i], VortexLineType::NewShed});
            }
        }
    }

} // namespace fvw