#include "wake.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

namespace fvw
{
    // Biot-Savart function
    // Optimize for blade structure
    void computeInducedVelocity(std::vector<Vec3> &inducedVelocities, const std::vector<Vec3> &targetPoints,
                                const Wake &wake, int timestep, const TurbineParams &turbineParams, double cutOff)
    {
        inducedVelocities.assign(targetPoints.size(), Vec3(0.0, 0.0, 0.0));

        if (timestep >= wake.bladeWakes.size() || wake.bladeWakes[timestep].empty())
        {
            std::cerr << "Warning: Trying to compute induced velocity for an invalid or empty timestep: " << timestep << std::endl;
            return; // 或者抛出异常
        }

        for (size_t p_idx = 0; p_idx < targetPoints.size(); ++p_idx)
        {
            Vec3 p = targetPoints[p_idx];
            Vec3 total_vel_at_p(0.0, 0.0, 0.0);

            // 遍历该时间步的所有叶片
            for (int b = 0; b < wake.nBlades; ++b)
            {
                const BladeWake &bladeWake = wake.getBladeWake(timestep, b);
                const auto &nodes = bladeWake.nodes; // 当前叶片的节点
                const auto &lines = bladeWake.lines; // 当前叶片的线段

                // 遍历当前叶片的所有线段
                for (const auto &line : lines)
                {
                    // 检查索引是否有效
                    if (line.startNodeIdx >= nodes.size() || line.endNodeIdx >= nodes.size() || line.startNodeIdx < 0 || line.endNodeIdx < 0)
                    {
                        std::cerr << "Warning: Invalid node index encountered in computeInducedVelocity for blade " << b << ", timestep " << timestep << std::endl;
                        continue; // 跳过这条无效线段
                    }

                    Vec3 x1 = nodes[line.startNodeIdx].position;
                    Vec3 x2 = nodes[line.endNodeIdx].position;
                    double gamma = line.gamma;

                    // --- 执行 Biot-Savart 计算 ---

                    Vec3 l_vec = x2 - x1;
                    double l_squared = l_vec.x * l_vec.x + l_vec.y * l_vec.y + l_vec.z * l_vec.z;

                    if (l_squared < 1e-12)
                    {
                        std::cerr << "Warning: The line is smaller than 1e-12" << std::endl;
                        continue; // 跳过长度为零的线段
                    }

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

                    // --- 结束 Biot-Savart 计算 ---
                    total_vel_at_p = total_vel_at_p + vel_contrib;
                }
            }

            inducedVelocities[p_idx] = total_vel_at_p;
        }
    }

    // --- initializeWake ---
    // The first wake
    // 初始化 t=0 时刻的尾迹结构 (附着涡 + 第一层脱落涡)
    void InitializeWakeStructure(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                                 const TurbineParams &turbineParams, const PositionData &pos, double dt)
    {
        bool if_verbose = false;

        // 清理并确保 t=0 和 t=1 的结构存在
        wake.bladeWakes.clear();
        wake.ensureTimeStepExists(0); // 创建 bladeWakes[0]
        wake.ensureTimeStepExists(1); // 创建 bladeWakes[1] 的空结构

        Vec3 initialFreeStreamVel = Vec3(turbineParams.windSpeed, 0.0, 0.0); // 简化假设

        // --- 填充 t=0 时刻的尾迹 ---
        std::cout << "Initializing wake structure at t=0..." << std::endl;

        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &currentBladeWake = wake.getBladeWake(0, b); // 获取 t=0, 叶片 b 的数据结构

            // 1. 添加 t=0 的节点 (叶片上的点 + 初始尾迹点)
            // 通常，FVW模型在叶片上定义节点 (如1/4弦长点)，这些点也作为尾迹的起点
            for (int i = 0; i < wake.nTrail; ++i) // nTrail = nShed + 1
            {
                // 叶片上的节点 (例如，1/4弦长点，作为附着涡的节点)
                Vec3 boundPos = pos.quarterAt(b, 0, i);                                                            // 获取 t=0 时刻的位置
                currentBladeWake.boundNodeIndices[i] = currentBladeWake.addNode({boundPos, initialFreeStreamVel}); // 添加节点并获取索引

                // 初始尾迹节点 (例如，后缘点，紧随叶片之后)
                // 这里用 trailAt，假设它代表后缘位置
                Vec3 trailPos = pos.trailAt(b, 0, i);
                currentBladeWake.trailNodeIndices[i] = currentBladeWake.addNode({trailPos, initialFreeStreamVel}); // 添加节点并获取索引
            }

            // 2. 计算 t=0 的涡线强度 (基于初始 BEM 结果)
            std::vector<double> gamma_bound(wake.nShed);
            std::vector<double> gamma_trail(wake.nTrail);

            for (int i = 0; i < wake.nShed; ++i) // 附着涡段数
            {
                // 使用 Perf 中 t=0 的 Cl 值
                double cl = perf.clAt(b, 0, i);
                double chord = geom.chordShedding[i]; // 使用控制点/涡段对应的弦长
                gamma_bound[i] = 0.5 * turbineParams.windSpeed * chord * cl;
            }

            // 计算初始尾迹涡强度
            gamma_trail[0] = gamma_bound[0]; // 翼根处
            for (int i = 1; i < wake.nShed; ++i)
            {
                gamma_trail[i] = gamma_bound[i] - gamma_bound[i - 1];
            }
            gamma_trail[wake.nShed] = -gamma_bound[wake.nShed - 1]; // 翼尖处

            // 3. 添加 t=0 的涡线
            // 添加附着涡线 (Bound)
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.boundNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i + 1];
                currentBladeWake.addLine({startIdx, endIdx, gamma_bound[i], VortexLineType::Bound});
            }

            // 添加初始尾迹涡线 (Trailing) - 连接叶片和第一层尾迹节点
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i];
                currentBladeWake.addLine({startIdx, endIdx, gamma_trail[i], VortexLineType::Trailing});
            }

            // 添加初始脱落涡线 (Shed) - 连接相邻的初始尾迹节点
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.trailNodeIndices[i + 1];
                currentBladeWake.addLine({startIdx, endIdx, -gamma_bound[i], VortexLineType::Shed});
            }
        }

        if (if_verbose)
        {
            std::cout << "Wake t=0 initialized for " << wake.nBlades << " blades." << std::endl;
            // 可以添加更详细的输出，例如第一个叶片的节点和线段数量
            if (wake.nBlades > 0)
            {
                std::cout << "  Blade 0 (t=0): " << wake.getNodes(0, 0).size() << " nodes, "
                          << wake.getLines(0, 0).size() << " lines." << std::endl;
            }
        }

        // --- 计算 t=0 节点的初始速度 (诱导速度 + 自由流) ---
        std::cout << "Computing initial velocities for t=0 nodes..." << std::endl;
        UpdateWakeVelocities(wake, turbineParams, 0);
    }

    // ------------------------

    void UpdateWakeVelocities(Wake &wake, const TurbineParams &turbineParams, int timestep)
    {
        // 计算诱导速度，叠加到节点速度     
        Vec3 initialFreeStreamVel = Vec3(turbineParams.windSpeed, 0.0, 0.0);

        // 1. 收集所有节点的位置作为目标点
        std::vector<Vec3> allNodePositions;
        std::vector<std::pair<int, int>> nodeGlobalIndices; // 存储 (blade, local_node_idx)
        for (int b = 0; b < wake.nBlades; ++b)
        {
            const auto &nodes = wake.getNodes(timestep, b);
            for (size_t n_idx = 0; n_idx < nodes.size(); ++n_idx)
            {
                allNodePositions.push_back(nodes[n_idx].position);
                nodeGlobalIndices.push_back({b, static_cast<int>(n_idx)});
            }
        }

        if (allNodePositions.empty())
        {
            std::cout << "No nodes found at timestep " << timestep << " to update velocities." << std::endl;
            return;
        }

        // 2. 计算这些目标点的诱导速度
        std::vector<Vec3> inducedVelocities;
        computeInducedVelocity(inducedVelocities, allNodePositions, wake, timestep, turbineParams);

        // 3. 更新 Wake 结构中对应节点的 velocity
        if (inducedVelocities.size() != allNodePositions.size())
        {
            std::cerr << "Error: Mismatch between number of target points and calculated induced velocities." << std::endl;
            return; // 或者抛出异常
        }

        for (size_t i = 0; i < inducedVelocities.size(); ++i)
        {
            int bladeIdx = nodeGlobalIndices[i].first;
            int localNodeIdx = nodeGlobalIndices[i].second;
            // 总速度 = 自由流速度 + 诱导速度
            wake.getNodes(timestep, bladeIdx)[localNodeIdx].velocity = initialFreeStreamVel + inducedVelocities[i];
        }
        std::cout << "Node velocities updated for timestep " << timestep << "." << std::endl;
    }

    // 把wake向前convect一步
    // 
    void AdvanceWakeStructure(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                              const TurbineParams &turbineParams, const PositionData &pos,
                              double dt, int currentTimestep)
    {
        // 确保 t=n 存在
        wake.ensureTimeStepExists(currentTimestep);

        Vec3 initialFreeStreamVel = Vec3(turbineParams.windSpeed, 0.0, 0.0);

        // 更新每个叶片的尾迹
        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);
            const BladeWake &prevBladeWake = wake.getBladeWake(currentTimestep - 1, b);

            // 0. 存储 t=n-1 的 Bound 涡量强度
            // Question：这一步会不会特别耗时？可能需要优化
            currentBladeWake.prevGammaBound.assign(wake.nShed, 0.0);
            for (const auto &line : prevBladeWake.lines)
            {
                if (line.type == VortexLineType::Bound)
                {
                    for (int i = 0; i < wake.nShed; ++i)
                    {
                        if (line.startNodeIdx == prevBladeWake.boundNodeIndices[i] &&
                            line.endNodeIdx == prevBladeWake.boundNodeIndices[i + 1])
                        {
                            currentBladeWake.prevGammaBound[i] = line.gamma;
                            break;
                        }
                    }
                }
            }

            // 1. 对流 t=n-1 的所有节点到 t=n
            currentBladeWake.nodes.clear();
            currentBladeWake.lines.clear();
            currentBladeWake.boundNodeIndices.resize(wake.nTrail, -1);
            currentBladeWake.trailNodeIndices.resize(wake.nTrail, -1);

            // 对流节点：保持 t=n-1 的索引顺序
            // 对流所有涡点 (这个需要调整一下,是对流所有的涡点)
            std::vector<int> oldToNewIdx(prevBladeWake.nodes.size(), -1);
            for (size_t i = 0; i < prevBladeWake.nodes.size(); ++i)
            {
                const VortexNode &node = prevBladeWake.nodes[i];
                Vec3 newPos = node.position + node.velocity * dt;
                int newIdx = currentBladeWake.addNode({newPos, initialFreeStreamVel});
                oldToNewIdx[i] = newIdx;

                // 更新 trailNodeIndices
                for (int j = 0; j < wake.nTrail; ++j)
                {
                    // 原来的bound vortex变成现在的trail
                    if (i == prevBladeWake.boundNodeIndices[j])
                    {
                        currentBladeWake.trailNodeIndices[j] = newIdx;
                    }
                }
            }

            // 2. 添加 t=n 的新附着涡节点
            std::vector<int> newBoundNodeIndices(wake.nTrail);
            for (int i = 0; i < wake.nTrail; ++i)
            {
                Vec3 boundPos = pos.quarterAt(b, currentTimestep, i);
                newBoundNodeIndices[i] = currentBladeWake.addNode({boundPos, initialFreeStreamVel});
                currentBladeWake.boundNodeIndices[i] = newBoundNodeIndices[i];
            }

            // 3. 添加 t=n 的涡量线(拓扑结构) - Gamma将在Kutta迭代中确定
            // 复制 t=n-1 的 Trailing 和 Shed 涡量线，更新索引
            for (const auto &line : prevBladeWake.lines)
            {
                if (line.type == VortexLineType::Trailing || line.type == VortexLineType::Shed)
                {
                    int newStartIdx = oldToNewIdx[line.startNodeIdx];
                    int newEndIdx = oldToNewIdx[line.endNodeIdx];
                    if (newStartIdx < 0 || newEndIdx < 0)
                    {
                        std::cerr << "Invalid mapped index for line at t=" << currentTimestep << std::endl;
                        continue;
                    }
                    currentBladeWake.addLine({newStartIdx, newEndIdx, line.gamma, line.type});
                }
            }

            // 新 Bound 涡量线
            // 需要改：应该用上一步的bound update -- 已改
            std::vector<double> gamma_bound(wake.nShed);
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = newBoundNodeIndices[i];
                int endIdx = newBoundNodeIndices[i + 1];
                gamma_bound[i] = currentBladeWake.prevGammaBound[i];
                currentBladeWake.addLine({startIdx, endIdx, gamma_bound[i], VortexLineType::Bound});
            }

            // 新 Trailing 涡量线
            std::vector<double> gamma_trail(wake.nTrail);
            gamma_trail[0] = gamma_bound[0];
            for (int i = 1; i < wake.nShed; ++i)
            {
                gamma_trail[i] = gamma_bound[i] - gamma_bound[i - 1];
            }
            gamma_trail[wake.nShed] = -gamma_bound[wake.nShed - 1];
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = newBoundNodeIndices[i]; 
                currentBladeWake.addLine({startIdx, endIdx, gamma_trail[i], VortexLineType::Trailing});
            }

            // 新 Shed 涡量线
            // 需要改：应该是这一步bound和上一步Bound的差分
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.trailNodeIndices[i + 1];
                // 这里不太对，需要改 -- 暂时初始化gamma为0.0，后面在kutta再更新？
                currentBladeWake.addLine({startIdx, endIdx, 0.0, VortexLineType::Shed});
            }

            // 接下来就应该调用kutta循环来更新gamma？
        }
    }

    // --------------------------------
    // This is the main function of initializing the wake
    void InitializeWake(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                        const TurbineParams &turbineParams, const PositionData &pos, double dt)
    {
        // t=0
        InitializeWakeStructure(wake, geom, perf, turbineParams, pos, dt);

        // --- 初始化 t=1 ---
        std::cout << "Preparing structure for t=1..." << std::endl;
        AdvanceWakeStructure(wake, geom, perf, turbineParams, pos, dt, 1);
    }

} // namespace fvw

// std::pair<double, double> interpolateClCd(int airfoilIdx, double aoa, std::vector<AirfoilData> &airfoils)
// {
//     if (airfoilIdx < 0 || airfoilIdx >= static_cast<int>(airfoils.size()))
//     {
//         throw std::invalid_argument("Invalid airfoilIdx");
//     }

//     double cl = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cl, aoa);
//     double cd = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cd, aoa);
//     return {cl, cd};
// }

// // 这里的输入的lines应该是只有lifting line 也就是Bound
// void kuttaJoukowskiIteration(Wake &wake, PerformanceData &perf, const BladeGeometry &geom, NodeAxes &axes,
//                              const TurbineParams &turbineParams, const PositionData &pos, VelBCS &velBCS, std::vector<AirfoilData> &airfoils)
// {

//     int max_iter_Kutta = 2;
//     int current_t = 1;

//     double relaxation_factor = 0.3;

//     std::vector<fvw::VortexLine> linesLiftingGamma(wake.nBlades * wake.nShed);
//     for (size_t i = 0; i < wake.lines[1].size(); ++i)
//     {
//         if (wake.lines[1][i].type == fvw::VortexLineType::Bound)
//         {
//             linesLiftingGamma.push_back(wake.lines[1][i]);
//         }
//     }

//     // 提取出Pos['bound']的值
//     std::vector<VortexNode> liftingLinesNode(turbineParams.nBlades * turbineParams.nSegments);
//     for (int b = 0; b < turbineParams.nBlades; ++b)
//     {
//         for (int i = 0; i < turbineParams.nSegments; ++i)
//         {

//             liftingLinesNode.push_back({pos.boundAt(b, 0, i)});
//         }
//     }

//     int nShed = turbineParams.nSegments;

//     for (int iter = 0; iter < max_iter_Kutta; ++iter)
//     {
//         std::vector<Vec3> vel_uind_ll(liftingLinesNode.size(), Vec3(0.0, 0.0, 0.0));
//         // 对每一个lifting line上的node，计算所有Line的contribution
//         // 需要debug: This is not node, it the middle of the line (need to re-define)
//         // Nodes debug完了，但是wake.nodes[current_t], wake.lines[current_t]的构造似乎有问题，导致了段错误，需要进一步debug
//         computeInducedVelocity(vel_uind_ll, liftingLinesNode, wake.nodes[current_t], wake.lines[current_t], turbineParams);

//         // 坐标转化
//         std::vector<Vec3> vel_tot(liftingLinesNode.size());
//         std::vector<double> Vinf(liftingLinesNode.size());
//         std::vector<int> NFoil(liftingLinesNode.size());
//         std::vector<double> chord(liftingLinesNode.size());

//         int node_idx_counter = 0;

//         // 伪代码
//         // 对于lifting line计算出induced velocity后
//         // 需要把该Induced velocity从惯性坐标系转换到叶片坐标系
//         // （目前的问题：如果用bxnAt等更新，vel_uind_ll是三个叶片都放在一起的，而axes.bxnAt是分叶片和时间的，无法一一对应地调用）
//         // 目前尝试根据Idx调用
//         // 再在叶片坐标系的基础上加上vel_blade，得到叶片坐标系的总速度
//         for (int b = 0; b < turbineParams.nBlades; ++b)
//         {
//             for (int i = 0; i < nShed; ++i)
//             {
//                 Vec3 uind = vel_uind_ll[node_idx_counter];
//                 vel_tot[node_idx_counter].x = velBCS.at(b, current_t, i).x + axes.bxtAt(b, current_t, i).dot(uind);
//                 vel_tot[node_idx_counter].y = velBCS.at(b, current_t, i).y + axes.bytAt(b, current_t, i).dot(uind);
//                 vel_tot[node_idx_counter].z = velBCS.at(b, current_t, i).z + axes.bztAt(b, current_t, i).dot(uind);
//                 Vinf[node_idx_counter] = vel_tot[node_idx_counter].x * vel_tot[node_idx_counter].x + vel_tot[node_idx_counter].y * vel_tot[node_idx_counter].y + vel_tot[node_idx_counter].z * vel_tot[node_idx_counter].z;
//                 NFoil[node_idx_counter] = geom.airfoilIndex[i];
//                 chord[node_idx_counter] = geom.chordShedding[i];
//                 node_idx_counter++;
//             }
//         }

//         std::vector<double> aoa(liftingLinesNode.size());
//         std::vector<double> cl(liftingLinesNode.size());
//         std::vector<double> cd(liftingLinesNode.size());
//         std::vector<double> gamma_new(liftingLinesNode.size());
//         std::vector<double> dg(liftingLinesNode.size());
//         double max_dg = 0.0;

//         for (int i = 0; i < liftingLinesNode.size(); i++)
//         {
//             // 然后使用arctan计算出攻角
//             aoa[i] = std::atan2(-vel_tot[i].y, vel_tot[i].x) * 180.0 / M_PI;

//             // 在图表中插值提取出cl和cd
//             auto [cl_value, cd_value] = interpolateClCd(NFoil[i], aoa[i], airfoils);
//             cl[i] = cl_value;
//             cd[i] = cd_value;

//             // 计算bound vorticity
//             // 公式: 0.5*Vinf*chord*cl
//             gamma_new[i] = 0.5 * Vinf[i] * chord[i] * cl[i];
//             dg[i] = gamma_new[i] - linesLiftingGamma[i].gamma;

//             // 使用relaxation_factor更新bound vortex
//             linesLiftingGamma[i].gamma = linesLiftingGamma[i].gamma + relaxation_factor * dg[i];
//             max_dg = std::max(max_dg, std::abs(dg[i]));
//         }

//         // 更新trailing vortex和shed vortex，这里又涉及到一个问题就是，trail和shed要针对不同的blade计算。
//         // 我觉得为了更好分析，是不是还是应该分成不同blade的？

//         // 判断bound vortex的error是否小于tol_kutta,如果是就可以跳出循环

//         // 预期返回值: bound, trail, shed vorticity 和 perf中cl, cd, aoa
//     }
// }

// } // namespace fvw