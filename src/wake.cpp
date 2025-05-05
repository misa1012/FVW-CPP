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

        // 用于存储 t=0 时刻每个叶片上附着涡节点和尾迹节点的内部索引
        std::vector<std::vector<int>> boundNodeIndices_t0(wake.nBlades, std::vector<int>(wake.nTrail));
        std::vector<std::vector<int>> trailNodeIndices_t0(wake.nBlades, std::vector<int>(wake.nTrail));
        // 用于存储 t=0 时刻计算出的附着涡和尾迹涡强度
        std::vector<std::vector<double>> gamma_bound_t0(wake.nBlades, std::vector<double>(wake.nShed));
        std::vector<std::vector<double>> gamma_trail_t0(wake.nBlades, std::vector<double>(wake.nTrail));

        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &currentBladeWake_t0 = wake.getBladeWake(0, b); // 获取 t=0, 叶片 b 的数据结构

            // 1. 添加 t=0 的节点 (叶片上的点 + 初始尾迹点)
            // 通常，FVW模型在叶片上定义节点 (如1/4弦长点)，这些点也作为尾迹的起点
            for (int i = 0; i < wake.nTrail; ++i) // nTrail = nShed + 1
            {
                // 叶片上的节点 (例如，1/4弦长点，作为附着涡的节点)
                Vec3 boundPos = pos.quarterAt(b, 0, i);                                                    // 获取 t=0 时刻的位置
                boundNodeIndices_t0[b][i] = currentBladeWake_t0.addNode({boundPos, initialFreeStreamVel}); // 添加节点并获取索引

                // 初始尾迹节点 (例如，后缘点，紧随叶片之后)
                // 简单处理：初始尾迹点与叶片后缘重合，或稍微向下游移动一点点
                // 这里用 trailAt，假设它代表后缘位置
                Vec3 trailPos = pos.trailAt(b, 0, i);
                // 或者稍微向下游移动: trailPos = pos.trailAt(b, 0, i) + initialFreeStreamVel * dt * 0.01; // 极小位移
                trailNodeIndices_t0[b][i] = currentBladeWake_t0.addNode({trailPos, initialFreeStreamVel}); // 添加节点并获取索引
            }

            // 2. 计算 t=0 的涡线强度 (基于初始 BEM 结果)
            for (int i = 0; i < wake.nShed; ++i) // 附着涡段数
            {
                // 使用 Perf 中 t=0 的 Cl 值
                double cl = perf.clAt(b, 0, i);
                double chord = geom.chordShedding[i]; // 使用控制点/涡段对应的弦长
                gamma_bound_t0[b][i] = 0.5 * turbineParams.windSpeed * chord * cl;
            }

            // 计算初始尾迹涡强度
            gamma_trail_t0[b][0] = gamma_bound_t0[b][0]; // 翼根处
            for (int i = 1; i < wake.nShed; ++i)
            {
                gamma_trail_t0[b][i] = gamma_bound_t0[b][i] - gamma_bound_t0[b][i - 1];
            }
            gamma_trail_t0[b][wake.nShed] = -gamma_bound_t0[b][wake.nShed - 1]; // 翼尖处

            // 3. 添加 t=0 的涡线
            // 添加附着涡线 (Bound)
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = boundNodeIndices_t0[b][i];
                int endIdx = boundNodeIndices_t0[b][i + 1];
                currentBladeWake_t0.addLine({startIdx, endIdx, gamma_bound_t0[b][i], VortexLineType::Bound});
            }

            // 添加初始尾迹涡线 (Trailing) - 连接叶片和第一层尾迹节点
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = trailNodeIndices_t0[b][i]; // 从尾迹点开始
                int endIdx = boundNodeIndices_t0[b][i];   // 指向叶片点
                currentBladeWake_t0.addLine({startIdx, endIdx, gamma_trail_t0[b][i], VortexLineType::Trailing});
            }

            // 添加初始脱落涡线 (Shed) - 连接相邻的初始尾迹节点
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = trailNodeIndices_t0[b][i];
                int endIdx = trailNodeIndices_t0[b][i + 1];
                // Shed 涡的强度通常与附着涡相同，方向相反 (保证环量守恒)
                currentBladeWake_t0.addLine({startIdx, endIdx, -gamma_bound_t0[b][i], VortexLineType::Shed});
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

        InitializeWakeVelocities(wake, turbineParams);

        // --- 准备 t=1 的结构 ---
        std::cout << "Preparing structure for t=1..." << std::endl;
        // 用于存储 t=1 时刻叶片上新节点的内部索引
        std::vector<std::vector<int>> boundNodeIndices_t1(wake.nBlades, std::vector<int>(wake.nTrail));

        // 用于存储 t=0 节点对流到 t=1 后的新索引
        std::vector<std::vector<int>> convectedBoundNodeIndices_t1(wake.nBlades, std::vector<int>(wake.nTrail));
        std::vector<std::vector<int>> convectedTrailNodeIndices_t1(wake.nBlades, std::vector<int>(wake.nTrail));

        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &currentBladeWake_t1 = wake.getBladeWake(1, b);       // 获取 t=1, 叶片 b 的数据结构
            const BladeWake &currentBladeWake_t0 = wake.getBladeWake(0, b); // 获取 t=0 的数据用于对流

            // 1. 对流 t=0 的节点到 t=1
            for (int i = 0; i < wake.nTrail; ++i)
            {
                // 对流附着涡节点 (虽然它在叶片上，但概念上它也产生尾迹)
                const VortexNode &node_b0 = currentBladeWake_t0.nodes[boundNodeIndices_t0[b][i]];
                Vec3 newPos_b = node_b0.position + node_b0.velocity * dt;
                // 添加到 t=1 的节点列表，速度暂时设为自由流 (之后会被重新计算)
                convectedBoundNodeIndices_t1[b][i] = currentBladeWake_t1.addNode({newPos_b, initialFreeStreamVel});

                // 对流尾迹节点
                const VortexNode &node_t0 = currentBladeWake_t0.nodes[trailNodeIndices_t0[b][i]];
                Vec3 newPos_t = node_t0.position + node_t0.velocity * dt;
                convectedTrailNodeIndices_t1[b][i] = currentBladeWake_t1.addNode({newPos_t, initialFreeStreamVel});
            }

            // 2. 添加 t=1 时刻叶片上的新节点 (附着涡节点)
            for (int i = 0; i < wake.nTrail; ++i) {
                Vec3 boundPos_t1 = pos.quarterAt(b, 1, i); // 获取 t=1 时刻的位置
                boundNodeIndices_t1[b][i] = currentBladeWake_t1.addNode({boundPos_t1, initialFreeStreamVel});
           }

            // 3. 添加 t=1 的涡线 (拓扑结构) - Gamma将在Kutta迭代中确定
            //    此时 Gamma 可以暂时设为 0 或 t=0 的值作为初始猜测

            
        }
    }

    // ------------------------

    void InitializeWakeVelocities(Wake &wake, const TurbineParams &turbineParams)
    {
        // 计算诱导速度，叠加到节点速度
        // --- 计算 t=0 节点的初始速度 (诱导速度 + 自由流) ---
        std::cout << "Computing initial velocities for t=0 nodes..." << std::endl;

        int timestep = 0;
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


    // --------------------------------
    // This is the main function of initializing the wake
    void InitializeWake(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                        const TurbineParams &turbineParams, const PositionData &pos, double dt)
    {
        // t=0
        InitializeWakeStructure(wake, geom, perf, turbineParams, pos, dt);
    }

} // namespace fvw

//     // Update velocity
//     std::vector<Vec3> inducedVel(nodes.size());
//     std::vector<Vec3> newTrailPos(nodes.size());

//     computeInducedVelocity(inducedVel, wake.nodes[0], wake.nodes[0], wake.lines[0], turbineParams);

//     for (int i = 0; i < NumNode0; ++i)
//     {
//         nodes[i].velocity = nodes[i].velocity + inducedVel[i];
//         newTrailPos[i] = nodes[i].position + nodes[i].velocity * dt;
//     }

//     // if (if_verbose)
//     // {
//     //     int blade = 0;
//     //     int timestep = 0;
//     //     // 输出指定叶片和时间步的诱导速度 x 坐标
//     //     std::cout << std::fixed << std::setprecision(6);
//     //     std::cout << "[Induced Velocity] Blade " << blade << ", t=" << timestep
//     //               << ", inducedVel.x:" << std::endl;
//     //     for (int i = 0; i < wake.nTrail; ++i)
//     //     {
//     //         int controlIdx = controlNodeIdx[blade][i];
//     //         int trailIdx = trailNodeIdx[blade][i];
//     //         std::cout << "i=" << std::setw(2) << i
//     //                   << ", controlNodeIdx=" << std::setw(3) << controlIdx
//     //                   << ", inducedVel.x=" << std::setw(10) << inducedVel[controlIdx].x << std::endl;
//     //         std::cout << "i=" << std::setw(2) << i
//     //                   << ", trailNodeIdx=" << std::setw(3) << trailIdx
//     //                   << ", inducedVel.x=" << std::setw(10) << inducedVel[trailIdx].x << std::endl;
//     //     }

//     //     // 输出指定叶片的第一个节点的 velocity
//     //     int firstNodeIdx = controlNodeIdx[blade][0];
//     //     std::cout << "First node vel (Blade " << blade << ", t=" << timestep
//     //               << ")=" << to_string(nodes[firstNodeIdx].velocity) << std::endl;
//     // }

//     // 第一步对流：使用前向欧拉法更新 trailing vortex 节点
//     std::cout << "Advancing timestep = 1" << std::endl;

//     wake.nodes.emplace_back(); // 创建 t=1 的节点数组，这之后不能再引用lines了，因为很可能对wake已经重新分配了内存
//     std::vector<VortexNode> &nodesNext = wake.nodes[1];
//     int NumNode1 = wake.nBlades * 3 * wake.nTrail;
//     nodesNext.reserve(NumNode1);

//     for (int i = 0; i < NumNode0; ++i)
//     {
//         nodesNext.push_back({newTrailPos[i], initialVel, wake.nodes[0][i].idx});
//     }

//     // 对于之前的Node：
//     // Convect all nodes forward in time using Forward Euler
//     // 对流所有 t=0 节点（绑定和尾迹节点）并添加新绑定节点

//     std::vector<std::vector<int>> newBoundNodeIdx(wake.nBlades, std::vector<int>(wake.nTrail));
//     for (int b = 0; b < wake.nBlades; ++b)
//     {
//         for (int i = 0; i < wake.nTrail; ++i)
//         {
//             Vec3 newBoundPosT1 = pos.quarterAt(b, 1, i); // t=1 的 1/4 弦长位置
//             newBoundNodeIdx[b][i] = node_idx_counter;
//             nodesNext.push_back({newBoundPosT1, initialVel, node_idx_counter++});
//         }
//     }

//     // Update line
//     wake.lines.emplace_back();
//     std::vector<VortexLine> &linesNext = wake.lines[1];
//     int NumLine1 = wake.nBlades * (3 * wake.nShed + 2 * wake.nTrail);
//     linesNext.reserve(NumLine1);

//     for (int n = 0; n < NumLine0; ++n)
//     {
//         if (wake.lines[0][n].type != VortexLineType::Bound)
//         {
//             linesNext.push_back(wake.lines[0][n]);
//         }
//     }

//     // Update bound vortex
//     for (int b = 0; b < wake.nBlades; ++b)
//     {
//         for (int i = 0; i < wake.nShed; ++i)
//         {
//             int startIdx = newBoundNodeIdx[b][i];
//             int endIdx = newBoundNodeIdx[b][i + 1];
//             linesNext.push_back({startIdx, endIdx, gamma_bound[b][i], VortexLineType::Bound});
//         }
//     }

//     // Update trail vortex
//     for (int b = 0; b < wake.nBlades; ++b)
//     {
//         for (int i = 0; i < wake.nTrail; ++i)
//         {
//             int startIdx = boundNodeIdx[b][i];
//             int endIdx = newBoundNodeIdx[b][i];
//             linesNext.push_back({startIdx, endIdx, gamma_trail[b][i], VortexLineType::Trailing});
//         }
//     }

//     // Update shed vortex
//     for (int b = 0; b < wake.nBlades; ++b)
//     {
//         for (int i = 0; i < wake.nTrail; ++i)
//         {
//             int startIdx = boundNodeIdx[b][i];
//             int endIdx = boundNodeIdx[b][i + 1];
//             linesNext.push_back({startIdx, endIdx, -1 * gamma_bound[b][i], VortexLineType::NewShed});
//         }
//     }
// }

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