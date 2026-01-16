#include "wake.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <omp.h>
#include <algorithm>

namespace fvw
{
    // Biot-Savart function
    // Optimize for blade structure
    void computeInducedVelocity(std::vector<Vec3> &inducedVelocities, const std::vector<Vec3> &targetPoints,
                                const Wake &wake, int timestep, const TurbineParams &turbineParams, const BladeGeometry &geom, const SimParams &simParams)
    {
        bool enableOpenMP =
#ifdef NDEBUG
            true;
#else
            false;
#endif

        const size_t numTargetPoints = targetPoints.size(); // 缓存大小
        inducedVelocities.assign(numTargetPoints, Vec3(0.0, 0.0, 0.0));

        if (timestep >= wake.bladeWakes.size() || wake.bladeWakes[timestep].empty())
        {
            std::cerr << "Warning: Trying to compute induced velocity for an invalid or empty timestep: " << timestep << std::endl;
            return; // 或者抛出异常
        }

        double four_pi = 4.0 * M_PI;

#pragma omp parallel for if (enableOpenMP)
        for (size_t p_idx = 0; p_idx < numTargetPoints; ++p_idx)
        {
            const Vec3 &p = targetPoints[p_idx]; // 使用引用
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

                    const Vec3 &x1 = nodes[line.startNodeIdx].position; // 使用引用
                    const Vec3 &x2 = nodes[line.endNodeIdx].position;
                    double gamma = line.gamma;

                    // Vortex core model
                    double smoothing_term;
                    switch (simParams.coreType)
                    {
                    case VortexCoreType::VanGarrel:
                    {
                        Vec3 l_vec = x2 - x1;
                        double l_squared = l_vec.norm_squared();
                        if (l_squared < 1e-12)
                        {
                            std::cerr << "Warning: The line is smaller than 1e-12" << std::endl;
                            continue; // 跳过长度为零的线段
                        }
                        smoothing_term = simParams.cutoffParam * simParams.cutoffParam * l_squared;
                        break;
                    }
                    case VortexCoreType::ChordBasedCore:
                    {
                        // 跟local chord挂钩
                        int segment_idx = line.segment_index;
                        double chord_local = 0.0;
                        if (line.type == VortexLineType::Bound || line.type == VortexLineType::Shed)
                        {
                            // 对于附着涡(Bound)和脱落涡(Shed)，它们与叶片上的“叶素分段”直接相关。
                            // 因此，使用在分段中心定义的 chordShedding 是物理上最合理的。
                            if (segment_idx >= 0 && segment_idx < geom.chordShedding.size())
                            {
                                chord_local = geom.chordShedding[segment_idx];
                            }
                            else
                            {
                                // 如果索引无效（例如，对于从旧代码继承的、没有索引的远场shed涡），使用一个平均值作为后备。
                                chord_local = 3.0; // 假设的平均弦长
                            }
                        }
                        else if (line.type == VortexLineType::Trailing)
                        {
                            // 对于尾迹涡(Trailing)，它们是从叶片后缘的“节点”上脱落的。
                            // 因此，使用在后缘节点上定义的 chordTrailing 更为精确。
                            // 在你的结构中，trailing line的segment_index对应于它所连接的后缘节点的索引。
                            if (segment_idx >= 0 && segment_idx < geom.chordTrailing.size())
                            {
                                chord_local = geom.chordTrailing[segment_idx];
                            }
                            else
                            {
                                // 后备方案
                                chord_local = 3.0;
                            }
                        }
                        else
                        {
                            // 为未知类型的涡线提供一个默认值
                            chord_local = 3.0;
                        }
                        smoothing_term = simParams.cutoffParam * simParams.cutoffParam * chord_local * chord_local;
                        // std::cout << "SegIdx: " << segment_idx
                        //           << ", Chord: " << chord_local
                        //           << ", SmoothingTerm: " << smoothing_term
                        //           << std::endl;
                        break;
                    }
                    }

                    double coeff = gamma / four_pi; // 使用预计算的 M_PI

                    Vec3 r1 = p - x1;
                    Vec3 r2 = p - x2;
                    double r1_norm = r1.norm();
                    double r2_norm = r2.norm();
                    double r1_r2 = r1_norm * r2_norm;
                    double dot_r1_r2 = r1.dot(r2);
                    Vec3 cross_r1_r2 = r1.cross(r2);

                    double denominator = r1_r2 * (r1_r2 + dot_r1_r2) + smoothing_term;

                    double factor = coeff * (r1_norm + r2_norm) / denominator;

                    total_vel_at_p.x += cross_r1_r2.x * factor;
                    total_vel_at_p.y += cross_r1_r2.y * factor;
                    total_vel_at_p.z += cross_r1_r2.z * factor;
                }
            }

            inducedVelocities[p_idx] = total_vel_at_p;
        }
    }

    // --- initializeWake ---
    // The first wake
    // 初始化 t=0 时刻的尾迹结构 (附着涡 + 第一层脱落涡)
    void InitializeWakeStructure(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                                 const TurbineParams &turbineParams, const PositionData &pos, const SimParams &simParams)
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
            currentBladeWake.boundLineIndices.resize(wake.nShed, -1);
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.boundNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i + 1];
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_bound[i], VortexLineType::Bound});
                VortexLine &newLine = currentBladeWake.lines.back();
                newLine.initial_gamma = newLine.gamma;
                newLine.in_far_wake = false;
                newLine.segment_index = i;
                currentBladeWake.boundLineIndices[i] = lineIdx;
            }

            // 添加初始尾迹涡线 (Trailing) - 连接叶片和第一层尾迹节点
            currentBladeWake.trailingLineIndices.resize(wake.nTrail, -1);
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i];
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_trail[i], VortexLineType::Trailing});
                VortexLine &newLine = currentBladeWake.lines.back();
                newLine.initial_gamma = newLine.gamma;
                newLine.in_far_wake = false;
                newLine.segment_index = i;
                currentBladeWake.trailingLineIndices[i] = lineIdx;
            }

            // 添加初始脱落涡线 (Shed) - 连接相邻的初始尾迹节点
            currentBladeWake.shedLineIndices.resize(wake.nShed, -1);
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.trailNodeIndices[i + 1];
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, -gamma_bound[i], VortexLineType::Shed});
                VortexLine &newLine = currentBladeWake.lines.back();
                newLine.initial_gamma = newLine.gamma;
                newLine.in_far_wake = false;
                newLine.segment_index = i;
                currentBladeWake.shedLineIndices[i] = lineIdx;
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
        UpdateWakeVelocities(wake, turbineParams, 0, geom, simParams);
    }

    // ------------------------
    // 该函数主要是在调用Biot-savart law，根据gamma计算每个点的速度
    void UpdateWakeVelocities(Wake &wake, const TurbineParams &turbineParams, int timestep, const BladeGeometry &geom, const SimParams &simParams)
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
        computeInducedVelocity(inducedVelocities, allNodePositions, wake, timestep, turbineParams, geom, simParams);

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
    // 该function的作用：
    // 1. 建立新的一步（t=n）的wake structure, how to archieve: wake.getBladeWake(currentTimestep, b)
    // 2. 读取上一步的bound vortex strength并存下来
    // 3. 处理新一步的节点：1) 把上一步的node根据速度convect；2)添加新的bound vortex node在quater chord line
    // 4. 处理新一步的vortex：1）把上一步的trail和shed复制到这一步（对应的start和end idx更新）；
    //                       2）生成新的bound和trail，需要add line附着到这上面（gamma暂时附一个值,后面需要通过Kutta循环更新)
    //                       3）上一步的bound变成这一步的shed,但是不需要上一步的信息了(直接丢弃),在这一步重新定义.(其vortex strength也需要在kutta循环中更新)
    void AdvanceWakeStructure(Wake &wake, const BladeGeometry &geom,
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

            int count = 0;
            for (const auto &line : prevBladeWake.lines)
            {
                if (line.type == VortexLineType::Bound)
                {
                    currentBladeWake.prevGammaBound[count] = line.gamma;
                    count++;
                }
            }

            // 1. 对流 t=n-1 的所有节点到 t=n
            currentBladeWake.nodes.clear();
            currentBladeWake.lines.clear();
            currentBladeWake.boundNodeIndices.assign(wake.nTrail, -1); // Resize and initialize with -1
            currentBladeWake.trailNodeIndices.assign(wake.nTrail, -1); // Resize and initialize with -1

            // Resize line index vectors
            currentBladeWake.boundLineIndices.assign(wake.nShed, -1);
            currentBladeWake.trailingLineIndices.assign(wake.nTrail, -1);
            currentBladeWake.shedLineIndices.assign(wake.nShed, -1);

            // 对流节点：保持 t=n-1 的索引顺序
            // 对流所有涡点 (这个需要调整一下,是对流所有的涡点)
            // std::vector<int> oldToNewIdx(prevBladeWake.nodes.size(), -1);
            for (size_t i = 0; i < prevBladeWake.nodes.size(); ++i)
            {
                const VortexNode &node = prevBladeWake.nodes[i];
                Vec3 newPos = node.position + node.velocity * dt;
                int _ = currentBladeWake.addNode({newPos, initialFreeStreamVel});
                // int newIdx = currentBladeWake.addNode({newPos, initialFreeStreamVel});
                // oldToNewIdx[i] = newIdx;
            }

            // 2. 添加 t=n 的新附着涡节点W
            // 并且更新这个时间步的trailNodeIndices和boundNodeIndices
            for (int i = 0; i < wake.nTrail; ++i)
            {
                currentBladeWake.trailNodeIndices[i] = prevBladeWake.boundNodeIndices[i];

                Vec3 boundPos = pos.quarterAt(b, currentTimestep, i);
                int newIdx = currentBladeWake.addNode({boundPos, initialFreeStreamVel});
                currentBladeWake.boundNodeIndices[i] = newIdx;
            }

            // -------------以上是对node的处理
            // -------------以下是对line的处理

            // 3. 添加 t=n 的涡量线(拓扑结构) - Gamma将在Kutta迭代中确定
            // 复制 t=n-1 的 Trailing 和 Shed 涡量线，更新索引
            for (const auto &line : prevBladeWake.lines)
            {
                if (line.type == VortexLineType::Trailing || line.type == VortexLineType::Shed)
                {
                    int lineIdx = currentBladeWake.addLine({line.startNodeIdx, line.endNodeIdx, line.gamma, line.type});
                    VortexLine &newLine = currentBladeWake.lines.back();
                    newLine.initial_gamma = line.initial_gamma; // 继承上一时间步的初始gamma
                    newLine.in_far_wake = line.in_far_wake;     // 继承上一时间步的远场状态
                    newLine.segment_index = line.segment_index;
                }
            }

            // 新 Bound 涡量线
            // 需要改：应该用上一步的bound update -- 已改
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.boundNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i + 1];
                std::vector<double> initial_gamma(wake.nShed);
                if (i < currentBladeWake.prevGammaBound.size())
                {
                    initial_gamma[i] = currentBladeWake.prevGammaBound[i];
                }
                else
                {
                    std::cerr << "Warning: prevGammaBound index out of bounds for blade " << b << ", segment " << i << " at t=" << currentTimestep << ". Initializing new bound gamma to 0." << std::endl;
                }

                // Add the new bound line and store its index
                int newLineIdx = currentBladeWake.addLine({startIdx, endIdx, initial_gamma[i], VortexLineType::Bound});
                if (i < currentBladeWake.boundLineIndices.size())
                {
                    VortexLine &newLine = currentBladeWake.lines.back();
                    newLine.initial_gamma = newLine.gamma;
                    newLine.in_far_wake = false;
                    newLine.segment_index = i;
                    currentBladeWake.boundLineIndices[i] = newLineIdx;
                }
                else
                {
                    std::cerr << "Error: boundLineIndices out of bounds when storing new bound line index for blade " << b << ", segment " << i << " at t=" << currentTimestep << std::endl;
                }
            }

            // 新 Trailing 涡量线
            std::vector<double> gamma_trail(wake.nTrail);
            gamma_trail[0] = currentBladeWake.prevGammaBound[0];
            for (int i = 1; i < wake.nShed; ++i)
            {
                gamma_trail[i] = currentBladeWake.prevGammaBound[i] - currentBladeWake.prevGammaBound[i - 1];
            }
            gamma_trail[wake.nShed] = -currentBladeWake.prevGammaBound[wake.nShed - 1];

            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i];
                int newLineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_trail[i], VortexLineType::Trailing});
                if (i < currentBladeWake.trailingLineIndices.size())
                {
                    VortexLine &newLine = currentBladeWake.lines.back();
                    newLine.initial_gamma = newLine.gamma;
                    newLine.in_far_wake = false;
                    newLine.segment_index = i;
                    currentBladeWake.trailingLineIndices[i] = newLineIdx;
                }
                else
                {
                    std::cerr << "Error: trailingLineIndices out of bounds when storing new trailing line index for blade " << b << ", segment " << i << " at t=" << currentTimestep << std::endl;
                }
            }

            // 新 Shed 涡量线
            // 需要改：应该是这一步bound和上一步Bound的差分
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.trailNodeIndices[i + 1];
                //  暂时初始化gamma，后面在kutta再更新
                int newLineIdx = currentBladeWake.addLine({startIdx, endIdx, 0.0, VortexLineType::Shed});
                if (i < currentBladeWake.shedLineIndices.size())
                {
                    VortexLine &newLine = currentBladeWake.lines.back();
                    // Shed涡的强度在Kutta循环中确定，所以初始gamma可以为0
                    newLine.initial_gamma = 0.0;
                    newLine.in_far_wake = false;
                    newLine.segment_index = i;
                    currentBladeWake.shedLineIndices[i] = newLineIdx;
                }
                else
                {
                    std::cerr << "Error: shedLineIndices out of bounds when storing new shed line index for blade " << b << ", segment " << i << " at t=" << currentTimestep << std::endl;
                }
            }

            // 接下来就应该调用kutta循环来更新gamma？
            // The Kutta-Joukowski iteration will be called after AdvanceWakeStructure
            // to update the gamma values for the newly created lines.
        }

        std::cout << "Wake structure advanced to timestep " << currentTimestep << "." << std::endl;
    }

    // Kutta循环 update vortex strength
    void kuttaJoukowskiIteration(Wake &wake, PerformanceData &perf, const BladeGeometry &geom, NodeAxes &axes,
                                 const TurbineParams &turbineParams, const PositionData &pos, VelBCS &velBCS,
                                 std::vector<AirfoilData> &airfoils, const SimParams &simParams)
    {
        const int maxIterations = 200;
        const double convergenceThreshold = 1e-4;
        const double relaxationFactor = 0.3;
        bool if_verbose = false;

        // 获取当前时间步（假设为最新的时间步）
        int currentTimestep = wake.bladeWakes.size() - 1;

        if (currentTimestep < 0)
        {
            std::cerr << "Error: No wake data available for Kutta-Joukowski iteration." << std::endl;
            return;
        }

        std::cout << "Starting Kutta-Joukowski iteration for timestep " << currentTimestep << "..." << std::endl;

        double max_dg = 0.0; // 用于检查收敛的最大 Gamma 变化量
        // 开始循环
        for (int iter = 0; iter < maxIterations; ++iter)
        {
            if (if_verbose)
            {
                std::cout << "  Kutta-Joukowski Iteration " << iter + 1 << "/" << maxIterations << std::endl;
            }

            max_dg = 0.0;

            // --- 1. 计算当前附着涡控制点 ---
            // 控制点通常取翼型的 1/4 弦长点 bound
            std::vector<Vec3> controlPointPositions;
            // 用于记录每个控制点对应哪个叶片和哪个叶片节段 (blade_idx, segment_idx)
            std::vector<std::pair<int, int>> controlPointMapping;

            for (int b = 0; b < wake.nBlades; ++b)
            {
                // 遍历当前叶片的所有附着涡线段对应的控制点 (nShed 个)
                for (int i = 0; i < wake.nShed; ++i)
                {
                    // 获取控制点位置 (1/4 弦长点)
                    try
                    {
                        Vec3 boundPos = pos.boundAt(b, currentTimestep, i);
                        controlPointPositions.push_back(boundPos);
                        // 记录控制点对应的叶片和节段索引
                        controlPointMapping.push_back({b, i});
                    }
                    catch (const std::out_of_range &e)
                    {
                        std::cerr << "Error accessing bound point for blade " << b << ", segment " << i << " at timestep " << currentTimestep << ": " << e.what() << std::endl;
                        continue; // 跳过这个无效的控制点
                    }
                }
            }

            // 如果没有控制点，则无法进行迭代
            if (controlPointPositions.empty())
            {
                std::cerr << "Warning: No control points found for Kutta iteration at timestep " << currentTimestep << "." << std::endl;
                break; // 退出 Kutta 循环
            }

            // --- 2. 计算这些控制点上的诱导速度 ---
            // 诱导速度来自整个尾迹结构 (所有时间步，所有叶片的所有涡线)
            std::vector<Vec3> inducedVelocities;
            // computeInducedVelocity 函数已经设计为计算整个 Wake 在目标点上的诱导速度
            computeInducedVelocity(inducedVelocities, controlPointPositions, wake, currentTimestep, turbineParams, geom, simParams);

            // --- 3. 计算有效速度，攻角，所需的附着涡强度，并更新附着涡 ---

            // 存储本迭代计算出的新的附着涡强度，用于后续更新 Trailing 和 Shed 涡线
            std::vector<std::vector<double>> updatedBoundGammas(wake.nBlades, std::vector<double>(wake.nShed));

            for (size_t k = 0; k < controlPointPositions.size(); ++k)
            {
                // 获取当前控制点对应的叶片和节段索引
                int b = controlPointMapping[k].first;
                int i = controlPointMapping[k].second; // 节段索引 (0 到 nShed-1)

                // 获取当前叶片在该时间步的 BladeWake 结构（非 const 版本，以便修改线段 gamma）
                BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);

                // 获取当前控制点上的诱导速度 (来自整个尾迹)
                Vec3 uind = inducedVelocities[k];

                // uind 是在惯性系下。需要将 uind 转换到叶片坐标系。
                const Vec3 &bxn = axes.bxnAt(b, currentTimestep, i); // 叶片切向 (弦向) 轴
                const Vec3 &byn = axes.bynAt(b, currentTimestep, i); // 垂直于弦，在旋转面内 (影响攻角)
                const Vec3 &bzn = axes.bznAt(b, currentTimestep, i); // 径向轴
                Vec3 uind_blade_frame;
                uind_blade_frame.x = bxn.dot(uind); // 沿弦方向的分量
                uind_blade_frame.y = byn.dot(uind); // 垂直于弦方向的分量 (影响攻角)
                uind_blade_frame.z = bzn.dot(uind); // 沿径向的分量

                // 计算有效速度 (在叶片坐标系下)
                Vec3 vel_eff_blade_frame = velBCS.at(b, currentTimestep, i) + uind_blade_frame;

                // 存储诱导速度
                // perf.setInducedVelocityICSAt(b, currentTimestep, i) = uind;
                // perf.setInducedVelocityAt(b, currentTimestep, i) = uind_blade_frame;
                perf.setRelativeVelocityAt(b, currentTimestep, i) = vel_eff_blade_frame;

                // 计算有效速度的模长
                double Vinf_eff_squared = vel_eff_blade_frame.x * vel_eff_blade_frame.x +
                                          vel_eff_blade_frame.y * vel_eff_blade_frame.y +
                                          vel_eff_blade_frame.z * vel_eff_blade_frame.z;
                double Vinf_eff = std::sqrt(Vinf_eff_squared);

                // 计算攻角 (AoA)
                double aoa_rad = std::atan2(-vel_eff_blade_frame.y, vel_eff_blade_frame.x);
                double aoa_deg = aoa_rad * 180.0 / M_PI;

                // 获取对应节段的翼型数据索引
                int airfoilIdx = geom.airfoilIndex[i];
                double cl_value = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cl, aoa_deg);
                double cd_value = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cd, aoa_deg);

                // 更新 PerformanceData
                perf.setAoaAt(b, currentTimestep, i) = aoa_deg;
                perf.setClAt(b, currentTimestep, i) = cl_value;
                perf.setCdAt(b, currentTimestep, i) = cd_value;

                // 计算维持当前 Cl 所需的附着涡强度 (根据 Kutta-Joukowski 定理)
                // 公式: Gamma = 0.5 * Rho * Vinf_eff * Chord * Cl
                double chord = geom.chordShedding[i];                      // 获取当前节段的弦长
                double gamma_required = 0.5 * Vinf_eff * chord * cl_value; // 考虑密度

                int boundLineIdx = currentBladeWake.boundLineIndices[i];
                VortexLine &boundLine = currentBladeWake.lines[boundLineIdx]; // 获取对线段的引用以便修改其 gamma

                // 使用松弛因子更新附着涡强度
                double gamma_current = boundLine.gamma;
                double dg = gamma_required - gamma_current;
                boundLine.gamma = gamma_current + relaxationFactor * dg;
                boundLine.initial_gamma = boundLine.gamma;

                // 存储到 PerformanceData
                perf.setBoundGammaAt(b, currentTimestep, i) = boundLine.gamma;

                // 存储本迭代计算出的新的 Bound Gamma，用于后续更新 Trailing/Shed
                if (i < updatedBoundGammas[b].size())
                {
                    updatedBoundGammas[b][i] = boundLine.gamma;
                }
                else
                {
                    std::cerr << "Error: updatedBoundGammas out of bounds for segment " << i << " on blade " << b << " at timestep " << currentTimestep << std::endl;
                }

                // 更新最大 Gamma 变化量，用于收敛判断
                max_dg = std::max(max_dg, std::abs(dg));
            } // 结束遍历控制点 (叶片节段)

            // --- 4. 根据更新后的附着涡强度，更新脱落涡 (Trailing) 和分离涡 (Shed) 的强度 ---
            // 注意：这里使用本迭代计算出的 updatedBoundGammas 来更新 Trailing 和 Shed 涡量。
            // Shed 涡量需要使用上一时间步的 Bound Gamma (prevGammaBound) 和当前时间步的 Bound Gamma。
            for (int b = 0; b < wake.nBlades; ++b)
            {
                BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);
                // Update Trailing 涡线 (nTrail 个)
                // Trailing line i 连接 trailNodeIndices[i] 到 boundNodeIndices[i]
                // 其强度等于其连接的两个附着点上附着涡的强度差。
                // Trailing 涡段 i 的强度 = Gamma_Bound[i] - Gamma_Bound[i-1] (对于 i > 0)
                // Trailing 涡段 0 (翼根) 的强度 = Gamma_Bound[0]
                // Trailing 涡段 nShed (翼尖) 的强度 = -Gamma_Bound[nShed-1]
                // 注意索引对应关系：nTrail = nShed + 1

                for (int i = 0; i < wake.nTrail; ++i) // 遍历 nTrail 个 Trailing 涡段 (0 到 nShed)
                {
                    // Find the current Trailing vortex line segment using its stored index
                    if (i < 0 || i >= currentBladeWake.trailingLineIndices.size())
                    {
                        std::cerr << "Error: trailingLineIndices out of bounds for segment " << i << " on blade " << b << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    int trailingLineIdx = currentBladeWake.trailingLineIndices[i];
                    if (trailingLineIdx == -1 || trailingLineIdx >= currentBladeWake.lines.size())
                    {
                        std::cerr << "Error: Invalid trailing line index " << trailingLineIdx << " in lines vector for blade " << b << ", segment " << i << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    VortexLine &trailingLine = currentBladeWake.lines[trailingLineIdx];

                    double gamma_trailing = 0.0; // Default value

                    if (i == 0) // 翼根处的 Trailing 涡段
                    {
                        // 强度等于第一个附着涡段的强度
                        if (wake.nShed > 0 && 0 < updatedBoundGammas[b].size())
                        {
                            gamma_trailing = updatedBoundGammas[b][0];
                        }
                        else if (wake.nShed == 0)
                        {
                            // No shed segments, no bound gamma to shed from
                            gamma_trailing = 0.0;
                        }
                        else
                        {
                            std::cerr << "Error: Index out of bounds accessing updatedBoundGammas[b][0] for blade " << b << ", segment " << i << ", timestep " << currentTimestep << "." << std::endl;
                        }
                    }
                    else if (i < wake.nShed) // 中间的 Trailing 涡段
                    {
                        // 强度等于其连接的两个附着涡段的强度差
                        if (i < updatedBoundGammas[b].size() && (i - 1) < updatedBoundGammas[b].size())
                        {
                            gamma_trailing = updatedBoundGammas[b][i] - updatedBoundGammas[b][i - 1];
                        }
                        else
                        {
                            std::cerr << "Error: Index out of bounds accessing updatedBoundGammas for blade " << b << ", segment " << i << ", timestep " << currentTimestep << "." << std::endl;
                        }
                    }
                    else if (i == wake.nShed) // 翼尖处的 Trailing 涡段
                    {
                        // 强度等于最后一个附着涡段强度的负值
                        if (wake.nShed > 0 && (wake.nShed - 1) < updatedBoundGammas[b].size())
                        {
                            gamma_trailing = -updatedBoundGammas[b][wake.nShed - 1];
                        }
                        else if (wake.nShed == 0)
                        {
                            // No shed segments
                            gamma_trailing = 0.0;
                        }
                        else
                        {
                            std::cerr << "Error: Index out of bounds accessing updatedBoundGammas for blade " << b << ", segment " << i << ", timestep " << currentTimestep << "." << std::endl;
                        }
                    }
                    else
                    {
                        std::cerr << "Error: Unexpected index i=" << i << " when updating trailing vortex for blade " << b << ", timestep " << currentTimestep << ". nTrail=" << wake.nTrail << std::endl;
                    }
                    trailingLine.gamma = gamma_trailing; // 更新 Trailing 涡线强度
                    trailingLine.initial_gamma = trailingLine.gamma;
                }

                // Update Shed 涡线 (nShed 个)
                // Shed line i 连接 trailNodeIndices[i] 到 trailNodeIndices[i+1]
                // 其强度等于上一时间步该位置对应的附着涡强度与当前时间步该位置对应的附着涡强度的差。
                // Shed 涡段 i 的强度 = Gamma_Bound[i] (t=n-1) - Gamma_Bound[i] (t=n)
                for (int i = 0; i < wake.nShed; ++i) // 遍历 nShed 个 Shed 涡段 (0 到 nShed-1)
                {
                    // Find the current Shed vortex line segment using its stored index
                    if (i < 0 || i >= currentBladeWake.shedLineIndices.size())
                    {
                        std::cerr << "Error: shedLineIndices out of bounds for segment " << i << " on blade " << b << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    int shedLineIdx = currentBladeWake.shedLineIndices[i];
                    if (shedLineIdx == -1 || shedLineIdx >= currentBladeWake.lines.size())
                    {
                        std::cerr << "Error: Invalid shed line index " << shedLineIdx << " in lines vector for blade " << b << ", segment " << i << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    VortexLine &shedLine = currentBladeWake.lines[shedLineIdx];

                    // Get previous timestep's bound gamma for this segment (stored in prevGammaBound)
                    double prev_gamma_bound = 0.0;
                    if (i < currentBladeWake.prevGammaBound.size())
                    {
                        prev_gamma_bound = currentBladeWake.prevGammaBound[i];
                    }
                    else
                    {
                        std::cerr << "Error: prevGammaBound out of bounds for blade " << b << ", segment " << i << " at timestep " << currentTimestep << ". Using 0 for previous gamma." << std::endl;
                    }

                    // Get current timestep's bound gamma for this segment (from updatedBoundGammas)
                    double current_gamma_bound = 0.0;
                    if (i < updatedBoundGammas[b].size())
                    {
                        current_gamma_bound = updatedBoundGammas[b][i];
                    }
                    else
                    {
                        std::cerr << "Error: updatedBoundGammas out of bounds for blade " << b << ", segment " << i << " at timestep " << currentTimestep << ". Using 0 for current gamma." << std::endl;
                    }

                    // Calculate Shed vortex line strength
                    shedLine.gamma = prev_gamma_bound - current_gamma_bound; // Note the sign convention
                    shedLine.initial_gamma = shedLine.gamma;

                } // 结束遍历 Shed 涡段
            }

            // --- 5. 检查收敛条件 ---
            if (max_dg < convergenceThreshold)
            {
                std::cout << "  Kutta-Joukowski iteration converged after " << iter + 1 << " iterations. Max |dGamma| = " << max_dg << std::endl;
                break; // 收敛，退出迭代循环
            }
        }

        // 如果迭代次数达到最大但未收敛，输出警告
        if (max_dg >= convergenceThreshold)
        {
            std::cerr << "Warning: Kutta-Joukowski iteration did not converge after " << maxIterations << " iterations. Max |dGamma| = " << max_dg << std::endl;
        }

        // Kutta 循环结束后，所有涡线的强度已经更新。
        // 在主循环中，下一步应该是根据这些新的涡线强度，再次计算所有节点的诱导速度，
        // 然后更新节点位置，进入下一个时间步。
    }

    void ApplyGammaDecayAndRemoval(Wake &wake, int timestep, const TurbineParams &turbineParams, const SimParams &simParams)
    {

        // --- 1. 定义模型参数 ---
        const double x_far_wake_D = 5.0; // 定义从下游5倍直径处开始衰减
        const double x_far_wake_m = x_far_wake_D * turbineParams.rTip * 2.0;

        const double alpha = 0.2;                  // 耗散率系数α (这是一个需要你调整的关键参数!)
        const double gamma_threshold_ratio = 0.01; // 当gamma衰减到初始值的1%时移除

        // --- 2. 遍历所有叶片和涡线，应用衰减 ---
        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &bladeWake = wake.getBladeWake(timestep, b);

            // a. 应用衰减
            for (VortexLine &line : bladeWake.lines)
            {
                // 附着涡不参与衰减
                if (line.type == VortexLineType::Bound)
                    continue;

                if (line.in_far_wake)
                {
                    line.gamma *= std::exp(-alpha * simParams.dt);
                }
                else
                {
                    const Vec3 &p1 = bladeWake.nodes[line.startNodeIdx].position;
                    const Vec3 &p2 = bladeWake.nodes[line.endNodeIdx].position;
                    double avg_x = (p1.x + p2.x) / 2.0;

                    if (avg_x > x_far_wake_m)
                    {
                        line.in_far_wake = true; // 标记它进入了远场
                    }
                }
            }

            // b. 移除弱涡
            auto &lines = bladeWake.lines;
            auto original_size = lines.size();

            // 使用C++ STL的 erase-remove idiom 高效地删除元素
            lines.erase(
                std::remove_if(lines.begin(), lines.end(),
                               [gamma_threshold_ratio](const VortexLine &line)
                               {
                                   if (!line.in_far_wake)
                                       return false;

                                   if (std::abs(line.initial_gamma) < 1e-9)
                                   {
                                       // 对于初始强度就接近0的涡（如某些shed涡），直接判断其绝对强度
                                       return std::abs(line.gamma) < 1e-7;
                                   }

                                   // 计算相对强度是否低于阈值
                                   return std::abs(line.gamma / line.initial_gamma) < gamma_threshold_ratio;
                               }),
                lines.end());

            if (original_size > lines.size())
            {
                std::cout << "  - Blade " << b << ": Removed " << original_size - lines.size() << " weak vortices from far-wake." << std::endl;
            }
        }
        // 注意：移除涡线后，与之关联的节点会变成“孤儿”节点。在这个模型中，我们可以暂时容忍
        // 这些孤儿节点的存在，因为它们不参与任何涡线的计算。一个更完整的实现需要垃圾回收机制。
    }

} // namespace fvw