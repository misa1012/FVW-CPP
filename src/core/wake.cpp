#include "core/wake.h"
#include "io/logger.h"
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
    // Optimized Biot-Savart function using Fast Filament approach
    void computeInducedVelocity(std::vector<Vec3> &inducedVelocities, const std::vector<Vec3> &targetPoints,
                                const Wake &wake, int timestep, const BladeGeometry &geom, const SimParams &simParams)
    {
        bool enableOpenMP =
#ifdef NDEBUG
            true;
#else
            false;
#endif

        const size_t numTargetPoints = targetPoints.size();
        inducedVelocities.assign(numTargetPoints, Vec3(0.0, 0.0, 0.0));

        if (timestep < 0 || static_cast<size_t>(timestep) >= wake.bladeWakes.size() || wake.bladeWakes[timestep].empty())
        {
            // std::cerr << "Warning: Trying to compute induced velocity for an invalid or empty timestep: " << timestep << std::endl;
            return;
        }

        // --- 1. Flatten / Linearize the Wake (Pre-computation) ---
        // This moves memory lookups and smoothing term calculations OUT of the N^2 loop.
        // It converts Array-of-Structures-with-Indirection to a flat Array-of-Structs.
        
        struct FastFilament
        {
            Vec3 x1;
            Vec3 x2;
            double gamma;
            double smoothing_term; 
        };

        // Estimate reserve size to avoid reallocation
        // Max possible = nBlades * (nShed + nTrail + nShed) ~ nBlades * 3 * nTrail
        size_t est_elements = wake.nBlades * (wake.nShed + wake.nShed + wake.nTrail);
        std::vector<FastFilament> fast_wake;
        fast_wake.reserve(est_elements);

        double four_pi = 4.0 * M_PI;
        double cutoff_sq = simParams.cutoffParam * simParams.cutoffParam;

        for (int b = 0; b < wake.nBlades; ++b)
        {
            const BladeWake &bladeWake = wake.getBladeWake(timestep, b);
            const auto &nodes = bladeWake.nodes;
            const auto &lines = bladeWake.lines;

            for (const auto &line : lines)
            {
                // Validity checks
                if (line.startNodeIdx < 0 || line.endNodeIdx < 0 ||
                    static_cast<size_t>(line.startNodeIdx) >= nodes.size() ||
                    static_cast<size_t>(line.endNodeIdx) >= nodes.size())
                {
                    continue;
                }

                // Skip weak vortices? (Optional optimization, strictness depends on needs)
                 if (std::abs(line.gamma) < 1e-9) continue; 

                const Vec3 &x1 = nodes[line.startNodeIdx].position;
                const Vec3 &x2 = nodes[line.endNodeIdx].position;
                
                // Pre-compute Geometry
                Vec3 l_vec = x2 - x1;
                double l_squared = l_vec.norm_squared();
                if (l_squared < 1e-12) continue; // Skip zero-length lines

                // Pre-compute Smoothing Term
                double smoothing = 0.0;
                if (simParams.coreType == VortexCoreType::VanGarrel)
                {
                    smoothing = cutoff_sq * l_squared;
                }
                else if (simParams.coreType == VortexCoreType::ChordBasedCore)
                {
                    int segment_idx = line.segment_index;
                    double chord_local = 3.0; // Default

                    if (line.type == VortexLineType::Bound || line.type == VortexLineType::Shed)
                    {
                        if (segment_idx >= 0 && static_cast<size_t>(segment_idx) < geom.chordShedding.size())
                            chord_local = geom.chordShedding[segment_idx];
                    }
                    else if (line.type == VortexLineType::Trailing)
                    {
                        if (segment_idx >= 0 && static_cast<size_t>(segment_idx) < geom.chordTrailing.size())
                            chord_local = geom.chordTrailing[segment_idx];
                    }
                    smoothing = cutoff_sq * chord_local * chord_local;
                }

                // Add to flat list
                fast_wake.push_back({x1, x2, line.gamma, smoothing});
            }
        }
        
        // --- 2. Parallel N-Body Interaction ---
        // Iterate over Target Points (Parallel)
        //    Iterate over FastFilaments (Serial, Linear Memory Access)

#pragma omp parallel for if (enableOpenMP) schedule(static)
        for (size_t p_idx = 0; p_idx < numTargetPoints; ++p_idx)
        {
            const Vec3 &p = targetPoints[p_idx];
            Vec3 total_vel_at_p(0.0, 0.0, 0.0);

            // Access raw pointer for potential compiler vectorization hints? 
            // Often iterators are optimized well enough in -O3.
            for (const auto &fil : fast_wake)
            {
                Vec3 r1 = p - fil.x1;
                Vec3 r2 = p - fil.x2;
                
                // Biot-Savart Kernel
                double r1_norm = r1.norm();
                double r2_norm = r2.norm();
                double r1_r2 = r1_norm * r2_norm;
                double dot_r1_r2 = r1.dot(r2);
                
                // Cross product
                // r1 x r2 = (p - x1) x (p - x2) = (p x p) - (p x x2) - (x1 x p) + (x1 x x2) 
                //         = 0 - p x x2 + p x x1 + x1 x x2
                //         = p x (x1 - x2) + x1 x x2
                // Actually r1 x r2 is standard.
                Vec3 cross_r1_r2 = r1.cross(r2);

                double denominator = r1_r2 * (r1_r2 + dot_r1_r2) + fil.smoothing_term;
                
                // Avoid divide by zero check? Smoothing term usually handles it, but robust code:
                if (denominator < 1e-12) continue;

                double coeff = fil.gamma / four_pi;
                double factor = coeff * (r1_norm + r2_norm) / denominator;

                total_vel_at_p.x += cross_r1_r2.x * factor;
                total_vel_at_p.y += cross_r1_r2.y * factor;
                total_vel_at_p.z += cross_r1_r2.z * factor;
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
        Logger::log("WAKE", "Initializing wake structure at t=0");

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
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_bound[i], VortexLineType::Bound, gamma_bound[i], false, i});
                currentBladeWake.lines.back();
                // newLine properties already set by initialization
                currentBladeWake.boundLineIndices[i] = lineIdx;
            }

            // 添加初始尾迹涡线 (Trailing) - 连接叶片和第一层尾迹节点
            currentBladeWake.trailingLineIndices.resize(wake.nTrail, -1);
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i];
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_trail[i], VortexLineType::Trailing, gamma_trail[i], false, i});
                currentBladeWake.lines.back();
                // properties set
                currentBladeWake.trailingLineIndices[i] = lineIdx;
            }

            // 添加初始脱落涡线 (Shed) - 连接相邻的初始尾迹节点
            currentBladeWake.shedLineIndices.resize(wake.nShed, -1);
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.trailNodeIndices[i + 1];
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, -gamma_bound[i], VortexLineType::Shed, -gamma_bound[i], false, i});
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
        Logger::log("VEL", "Computing initial velocities for t=0");
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
        computeInducedVelocity(inducedVelocities, allNodePositions, wake, timestep, geom, simParams);

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
        Logger::log("VEL", "Updated node velocities");
    }

    // 把wake向前convect一步
    // 该function的作用：
    // 1. 建立新的一步（t=n）的wake structure, how to archieve: wake.getBladeWake(currentTimestep, b)
    // 2. 读取上一步的bound vortex strength并存下来
    // 3. 处理新一步的节点：1) 把上一步的node根据速度convect；2)添加新的bound vortex node在quater chord line
    // 4. 处理新一步的vortex：1）把上一步的trail和shed复制到这一步（对应的start和end idx更新）；
    //                       2）生成新的bound和trail，需要add line附着到这上面（gamma暂时附一个值,后面需要通过Kutta循环更新)
    //                       3）上一步的bound变成这一步的shed,但是不需要上一步的信息了(直接丢弃),在这一步重新定义.(其vortex strength也需要在kutta循环中更新)
    void AdvanceWakeStructure(Wake &wake,
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
            for (size_t i = 0; i < prevBladeWake.nodes.size(); ++i)
            {
                const VortexNode &node = prevBladeWake.nodes[i];
                Vec3 newPos = node.position + node.velocity * dt;
                currentBladeWake.addNode({newPos, initialFreeStreamVel});
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

            // 3. 添加 t=n 的涡量线(拓扑结构) - Gamma将在Kutta迭代中确定
            // 复制 t=n-1 的 Trailing 和 Shed 涡量线，更新索引
            for (const auto &line : prevBladeWake.lines)
            {
                if (line.type == VortexLineType::Trailing || line.type == VortexLineType::Shed)
                {
                    currentBladeWake.addLine({line.startNodeIdx, line.endNodeIdx, line.gamma, line.type, line.initial_gamma, line.in_far_wake, line.segment_index});
                }
            }

            // 新 Bound 涡量线
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.boundNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i + 1];
                std::vector<double> initial_gamma(wake.nShed);
                if (static_cast<size_t>(i) < currentBladeWake.prevGammaBound.size())
                {
                    initial_gamma[i] = currentBladeWake.prevGammaBound[i];
                }
                else
                {
                    std::cerr << "Warning: prevGammaBound index out of bounds for blade " << b << ", segment " << i << " at t=" << currentTimestep << ". Initializing new bound gamma to 0." << std::endl;
                }

                // Add the new bound line and store its index
                int newLineIdx = currentBladeWake.addLine({startIdx, endIdx, initial_gamma[i], VortexLineType::Bound, initial_gamma[i], false, i});
                if (static_cast<size_t>(i) < currentBladeWake.boundLineIndices.size())
                {
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
                int newLineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_trail[i], VortexLineType::Trailing, gamma_trail[i], false, i});
                if (static_cast<size_t>(i) < currentBladeWake.trailingLineIndices.size())
                {
                    currentBladeWake.trailingLineIndices[i] = newLineIdx;
                }
                else
                {
                    std::cerr << "Error: trailingLineIndices out of bounds when storing new trailing line index for blade " << b << ", segment " << i << " at t=" << currentTimestep << std::endl;
                }
            }

            // 新 Shed 涡量线
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.trailNodeIndices[i + 1];
                //  暂时初始化gamma，后面在kutta再更新
                int newLineIdx = currentBladeWake.addLine({startIdx, endIdx, 0.0, VortexLineType::Shed, 0.0, false, i});
                if (static_cast<size_t>(i) < currentBladeWake.shedLineIndices.size())
                {
                    currentBladeWake.shedLineIndices[i] = newLineIdx;
                }
                else
                {
                    std::cerr << "Error: shedLineIndices out of bounds when storing new shed line index for blade " << b << ", segment " << i << " at t=" << currentTimestep << std::endl;
                }
            }
        }

        Logger::log("WAKE", "Advanced structure");
    }

    // Kutta循环 update vortex strength
    void kuttaJoukowskiIteration(Wake &wake, PerformanceData &perf, const BladeGeometry &geom, NodeAxes &axes,
                                 const PositionData &pos, VelBCS &velBCS,
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
            std::vector<Vec3> controlPointPositions;
            std::vector<std::pair<int, int>> controlPointMapping;

            for (int b = 0; b < wake.nBlades; ++b)
            {
                for (int i = 0; i < wake.nShed; ++i)
                {
                    try
                    {
                        Vec3 boundPos = pos.boundAt(b, currentTimestep, i);
                        controlPointPositions.push_back(boundPos);
                        controlPointMapping.push_back({b, i});
                    }
                    catch (const std::out_of_range &e)
                    {
                        std::cerr << "Error accessing bound point for blade " << b << ", segment " << i << " at timestep " << currentTimestep << ": " << e.what() << std::endl;
                        continue; 
                    }
                }
            }

            if (controlPointPositions.empty())
            {
                std::cerr << "Warning: No control points found for Kutta iteration at timestep " << currentTimestep << "." << std::endl;
                break; 
            }

            // --- 2. 计算这些控制点上的诱导速度 ---
            std::vector<Vec3> inducedVelocities;
            computeInducedVelocity(inducedVelocities, controlPointPositions, wake, currentTimestep, geom, simParams);

            // --- 3. 计算有效速度，攻角，所需的附着涡强度，并更新附着涡 ---
            std::vector<std::vector<double>> updatedBoundGammas(wake.nBlades, std::vector<double>(wake.nShed));

            for (size_t k = 0; k < controlPointPositions.size(); ++k)
            {
                int b = controlPointMapping[k].first;
                int i = controlPointMapping[k].second; 

                BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);

                Vec3 uind = inducedVelocities[k];
                const Vec3 &bxn = axes.bxnAt(b, currentTimestep, i); 
                const Vec3 &byn = axes.bynAt(b, currentTimestep, i); 
                const Vec3 &bzn = axes.bznAt(b, currentTimestep, i); 
                Vec3 uind_blade_frame;
                uind_blade_frame.x = bxn.dot(uind); 
                uind_blade_frame.y = byn.dot(uind); 
                uind_blade_frame.z = bzn.dot(uind); 

                Vec3 vel_eff_blade_frame = velBCS.at(b, currentTimestep, i) + uind_blade_frame;

                perf.setRelativeVelocityAt(b, currentTimestep, i) = vel_eff_blade_frame;

                double Vinf_eff_squared = vel_eff_blade_frame.x * vel_eff_blade_frame.x +
                                          vel_eff_blade_frame.y * vel_eff_blade_frame.y +
                                          vel_eff_blade_frame.z * vel_eff_blade_frame.z;
                double Vinf_eff = std::sqrt(Vinf_eff_squared);

                double aoa_rad = std::atan2(-vel_eff_blade_frame.y, vel_eff_blade_frame.x);
                double aoa_deg = aoa_rad * 180.0 / M_PI;

                int airfoilIdx = geom.airfoilIndex[i];
                double cl_value = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cl, aoa_deg);
                double cd_value = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cd, aoa_deg);

                perf.setAoaAt(b, currentTimestep, i) = aoa_deg;
                perf.setClAt(b, currentTimestep, i) = cl_value;
                perf.setCdAt(b, currentTimestep, i) = cd_value;

                double chord = geom.chordShedding[i];
                double gamma_required = 0.5 * Vinf_eff * chord * cl_value; 

                int boundLineIdx = currentBladeWake.boundLineIndices[i];
                VortexLine &boundLine = currentBladeWake.lines[boundLineIdx]; 

                double gamma_current = boundLine.gamma;
                double dg = gamma_required - gamma_current;
                boundLine.gamma = gamma_current + relaxationFactor * dg;
                boundLine.initial_gamma = boundLine.gamma;

                perf.setBoundGammaAt(b, currentTimestep, i) = boundLine.gamma;

                if (static_cast<std::size_t>(i) < updatedBoundGammas[b].size())
                {
                    updatedBoundGammas[b][i] = boundLine.gamma;
                }
                else
                {
                    std::cerr << "Error: updatedBoundGammas out of bounds for segment " << i << " on blade " << b << " at timestep " << currentTimestep << std::endl;
                }

                max_dg = std::max(max_dg, std::abs(dg));
            } 

            // --- 4. 更新脱落涡 (Trailing) 和分离涡 (Shed) ---
            for (int b = 0; b < wake.nBlades; ++b)
            {
                BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);

                for (int i = 0; i < wake.nTrail; ++i) 
                {
                    if (i < 0 || static_cast<size_t>(i) >= currentBladeWake.trailingLineIndices.size())
                    {
                        std::cerr << "Error: trailingLineIndices out of bounds for segment " << i << " on blade " << b << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    int trailingLineIdx = currentBladeWake.trailingLineIndices[i];
                    if (trailingLineIdx == -1 || static_cast<size_t>(trailingLineIdx) >= currentBladeWake.lines.size())
                    {
                        std::cerr << "Error: Invalid trailing line index " << trailingLineIdx << " in lines vector for blade " << b << ", segment " << i << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    VortexLine &trailingLine = currentBladeWake.lines[trailingLineIdx];

                    double gamma_trailing = 0.0; 

                    if (i == 0) 
                    {
                        if (wake.nShed > 0 && 0 < updatedBoundGammas[b].size())
                        {
                            gamma_trailing = updatedBoundGammas[b][0];
                        }
                    }
                    else if (i < wake.nShed) 
                    {
                        if (static_cast<size_t>(i) < updatedBoundGammas[b].size() && static_cast<size_t>(i - 1) < updatedBoundGammas[b].size())
                        {
                            gamma_trailing = updatedBoundGammas[b][i] - updatedBoundGammas[b][i - 1];
                        }
                    }
                    else if (i == wake.nShed) 
                    {
                        if (wake.nShed > 0 && static_cast<size_t>(wake.nShed - 1) < updatedBoundGammas[b].size())
                        {
                            gamma_trailing = -updatedBoundGammas[b][wake.nShed - 1];
                        }
                    }
                    trailingLine.gamma = gamma_trailing; 
                    trailingLine.initial_gamma = trailingLine.gamma;
                }

                for (int i = 0; i < wake.nShed; ++i) 
                {
                    if (i < 0 || static_cast<size_t>(i) >= currentBladeWake.shedLineIndices.size())
                    {
                        std::cerr << "Error: shedLineIndices out of bounds for segment " << i << " on blade " << b << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    int shedLineIdx = currentBladeWake.shedLineIndices[i];
                    if (shedLineIdx == -1 || static_cast<size_t>(shedLineIdx) >= currentBladeWake.lines.size())
                    {
                        std::cerr << "Error: Invalid shed line index " << shedLineIdx << " in lines vector for blade " << b << ", segment " << i << " at timestep " << currentTimestep << std::endl;
                        continue;
                    }
                    VortexLine &shedLine = currentBladeWake.lines[shedLineIdx];

                    double prev_gamma_bound = 0.0;
                    if (static_cast<size_t>(i) < currentBladeWake.prevGammaBound.size())
                    {
                        prev_gamma_bound = currentBladeWake.prevGammaBound[i];
                    }

                    double current_gamma_bound = 0.0;
                    if (static_cast<std::size_t>(i) < updatedBoundGammas[b].size())
                    {
                        current_gamma_bound = updatedBoundGammas[b][i];
                    }

                    shedLine.gamma = prev_gamma_bound - current_gamma_bound; 
                    shedLine.initial_gamma = shedLine.gamma;

                } 
            }

            if (max_dg < convergenceThreshold)
            {
                Logger::log("KUTTA", "Converged: " + std::to_string(iter + 1) + " iters, |dGamma| = " + std::to_string(max_dg));
                break; 
            }
        }

        if (max_dg >= convergenceThreshold)
        {
            Logger::log("KUTTA", "WARNING: Failed to converge after " + std::to_string(maxIterations) + 
                         " iters. Max |dGamma| = " + std::to_string(max_dg));
        }
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