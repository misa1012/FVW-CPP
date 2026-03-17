#include "core/wake.h"
#include "io/logger.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <omp.h>
#include <algorithm>
#include <thread>
#include <cstdlib>

namespace fvw
{
    namespace
    {
        // Shen-type single-factor tip correction (Wind Energy 2005; doi:10.1002/we.153)
        double shen_tip_factor(
            double r,
            double phi,
            const TurbineParams &turbineParams,
            double tip_speed_ratio_shen)
        {
            constexpr double pi = 3.14159265358979323846;
            const double sin_phi = std::max(std::abs(std::sin(phi)), 1e-6);
            const double r_safe = std::max(r, 1e-6);
            const double b = static_cast<double>(turbineParams.nBlades);

            // Same defaults as ALM branch for Shen-style correction.
            const double g = std::exp(-0.125 * ((b * tip_speed_ratio_shen) - 21.0)) + 0.1;
            const double ftip = (turbineParams.rTip - r_safe) / (r_safe * sin_phi + 1e-12);
            const double arg1 = -g * (b / 2.0) * ftip;
            const double arg2 = std::exp(std::max(std::min(30.0, arg1), -30.0));
            const double ftip_corr = (2.0 / pi) * std::acos(std::clamp(arg2, 0.0, 1.0));

            if (!std::isfinite(ftip_corr) || ftip_corr <= 0.0) return 1.0;
            return std::clamp(ftip_corr, 0.0, 1.0);
        }

        bool compute_phi(
            int b,
            int currentTimestep,
            int i,
            const PositionData &pos,
            const Vec3 &bxn,
            const Vec3 &byn,
            const Vec3 &bzn,
            const Vec3 &vel_eff_blade_frame,
            double &phi_out)
        {
            const Vec3 vel_eff_ics = bxn * vel_eff_blade_frame.x +
                                     byn * vel_eff_blade_frame.y +
                                     bzn * vel_eff_blade_frame.z;
            const Vec3 axis_ics(-1.0, 0.0, 0.0); // rotor thrust axis (opposite wind +x)
            Vec3 radial_ics = pos.boundAt(b, currentTimestep, i) - pos.hubAt(currentTimestep);
            radial_ics.x = 0.0; // project to rotor plane (y-z)
            const double radial_norm = radial_ics.norm();
            if (radial_norm <= 1e-12) return false;
            radial_ics = radial_ics * (1.0 / radial_norm);

            Vec3 tangential_ics = axis_ics.cross(radial_ics);
            const double tan_norm = tangential_ics.norm();
            if (tan_norm <= 1e-12) return false;
            tangential_ics = tangential_ics * (1.0 / tan_norm);

            const double v_axial = std::abs(vel_eff_ics.dot(axis_ics));
            const double v_tan = std::abs(vel_eff_ics.dot(tangential_ics));
            phi_out = std::atan2(v_axial, std::max(v_tan, 1e-12));
            return std::isfinite(phi_out);
        }

        double compute_tip_loss_factor(
            double r,
            double phi,
            const TurbineParams &turbineParams,
            TipLossModel model)
        {
            switch (model)
            {
                case TipLossModel::Off:
                    return 1.0;
                case TipLossModel::Shen:
                    return shen_tip_factor(r, phi, turbineParams, turbineParams.tsr);
            }
            return 1.0;
        }
    }

    // Optimized Biot-Savart function using Fast Filament approach
    void computeInducedVelocity(std::vector<Vec3> &inducedVelocities, const std::vector<Vec3> &targetPoints,
                                const Wake &wake, int timestep, const BladeGeometry &geom, const SimParams &simParams,
                                size_t startLineIdx, size_t endLineIdx)
    {
        bool enableOpenMP = true;
#ifndef _OPENMP
        enableOpenMP = false;
#endif
        
        static bool omp_status_printed = false;
        if (!omp_status_printed) {
            if (simParams.logVerbose) {
                if (enableOpenMP) {
                    // Auto-detect threads if environment variable is not set
                    if (std::getenv("OMP_NUM_THREADS") == nullptr) {
                        int hw_threads = std::thread::hardware_concurrency();
                        if (hw_threads > 0) {
                            omp_set_num_threads(hw_threads);
                            Logger::log("INFO", "OMP_NUM_THREADS not set. Auto-detected and set threads to: " + std::to_string(hw_threads));
                        }
                    }

                    int n_threads = omp_get_max_threads();
                    Logger::log("INFO", "OpenMP Enabled for Bio-Savart. Max threads: " + std::to_string(n_threads));
                } else {
                    Logger::log("WARN", "OpenMP Disabled (Compiler macro _OPENMP not found).");
                }

                // --- OpenMP Diagnostic Check ---
                // Verify if threads are actually spawning and running on distinct cores
                if (enableOpenMP) {
                    int max_t = omp_get_max_threads();
                    std::vector<int> thread_hits(max_t, 0);
                    #pragma omp parallel
                    {
                        #pragma omp atomic
                        thread_hits[omp_get_thread_num()]++;
                    }

                    int active_threads = 0;
                    for(int t : thread_hits) { if(t > 0) active_threads++; }
                    
                    Logger::log("INFO", "OpenMP Diagnostic: Requested " + std::to_string(max_t) + " threads, but " + std::to_string(active_threads) + " were active.");
                    if (active_threads < 2 && max_t > 1) {
                        Logger::log("WARN", "CRITICAL: OpenMP is active but only 1 thread is doing work! Check OMP_PROC_BIND or CPU affinity.");
                    }
                }
            }
            omp_status_printed = true;
        }

        const size_t numTargetPoints = targetPoints.size();
        
        // Always reset/overwrite output for the specified range
        inducedVelocities.assign(numTargetPoints, Vec3(0.0, 0.0, 0.0));

        if (timestep < 0 || static_cast<size_t>(timestep) >= wake.bladeWakes.size() || wake.bladeWakes[timestep].empty())
        {
            return;
        }

        // --- 1. Flatten / Linearize the Wake (Pre-computation) ---
        double t_prep_start = omp_get_wtime();
        // Precise reserve size calculation to avoid reallocations
        size_t est_elements = 0;
        for (int b = 0; b < wake.nBlades; ++b)
        {
             const auto &lines = wake.getBladeWake(timestep, b).lines;
             size_t limit = std::min(endLineIdx, lines.size());
             size_t start = std::min(startLineIdx, lines.size());
             if (limit > start) est_elements += (limit - start);
        }

        std::vector<double> x1x, x1y, x1z;
        std::vector<double> x2x, x2y, x2z;
        std::vector<double> gamma_over_4pi;
        std::vector<double> smoothing_term;
        x1x.reserve(est_elements);
        x1y.reserve(est_elements);
        x1z.reserve(est_elements);
        x2x.reserve(est_elements);
        x2y.reserve(est_elements);
        x2z.reserve(est_elements);
        gamma_over_4pi.reserve(est_elements);
        smoothing_term.reserve(est_elements);

        double cutoff_sq = simParams.cutoffParam * simParams.cutoffParam;
        double cutoff_4 = cutoff_sq * cutoff_sq;
        double four_pi = 4.0 * M_PI;

        for (int b = 0; b < wake.nBlades; ++b)
        {
            const BladeWake &bladeWake = wake.getBladeWake(timestep, b);
            const auto &nodes = bladeWake.nodes;
            const auto &lines = bladeWake.lines;

            size_t limit = std::min(endLineIdx, lines.size());
            size_t start = std::min(startLineIdx, lines.size());

            // Loop over specified range
            for (size_t i = start; i < limit; ++i)
            {
                const auto &line = lines[i];

                // Validity checks
                if (line.startNodeIdx < 0 || line.endNodeIdx < 0 ||
                    static_cast<size_t>(line.startNodeIdx) >= nodes.size() ||
                    static_cast<size_t>(line.endNodeIdx) >= nodes.size())
                {
                    continue;
                }

                if (std::abs(line.gamma) < 1e-9) continue; 

                const Vec3 &x1 = nodes[line.startNodeIdx].position;
                const Vec3 &x2 = nodes[line.endNodeIdx].position;
                
                Vec3 l_vec = x2 - x1;
                double l_squared = l_vec.norm_squared();
                if (l_squared < 1e-12) continue;

                // Pre-compute Smoothing Term
                double smoothing = 0.0;
                if (simParams.coreType == VortexCoreType::VanGarrel)
                {
                    // Van Garrel: epsilon = cutoffParam * L, smoothing = epsilon^2 (matches original formula form)
                    smoothing = cutoff_sq * l_squared;
                }
                else if (simParams.coreType == VortexCoreType::VanGarrelUnitConsistent)
                {
                    // Van Garrel unit-consistent: epsilon = cutoffParam * L, smoothing = epsilon^4
                    double l4 = l_squared * l_squared;
                    smoothing = cutoff_4 * l4;
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
                    // Chord-based: epsilon = cutoffParam * c, smoothing = epsilon^2 (matches original formula form)
                    double c2 = chord_local * chord_local;
                    smoothing = cutoff_sq * c2;
                }

                x1x.push_back(x1.x);
                x1y.push_back(x1.y);
                x1z.push_back(x1.z);
                x2x.push_back(x2.x);
                x2y.push_back(x2.y);
                x2z.push_back(x2.z);
                gamma_over_4pi.push_back(line.gamma / four_pi);
                smoothing_term.push_back(smoothing);
            }
        }
        double t_prep_end = omp_get_wtime();
        
        // --- 2. Parallel N-Body Interaction ---
        // HEURISTIC: Disable OpenMP for small workloads
        // Overhead of thread creation can dominate small loops.
        // Threshold ~100k interactions is safer.
        size_t workload = numTargetPoints * x1x.size();
        if (enableOpenMP && workload < 100000) {
             enableOpenMP = false;
        }

        double t_calc_start = omp_get_wtime();
#pragma omp parallel for if (enableOpenMP) schedule(static)
        for (size_t p_idx = 0; p_idx < numTargetPoints; ++p_idx)
        {
            const Vec3 &p = targetPoints[p_idx];
            Vec3 total_vel_at_p(0.0, 0.0, 0.0);

            const size_t n_fil = x1x.size();
            for (size_t f = 0; f < n_fil; ++f)
            {
                const Vec3 r1(p.x - x1x[f], p.y - x1y[f], p.z - x1z[f]);
                const Vec3 r2(p.x - x2x[f], p.y - x2y[f], p.z - x2z[f]);

                const double r1_sq = r1.norm_squared();
                const double r2_sq = r2.norm_squared();
                if (r1_sq < 1e-20 || r2_sq < 1e-20) continue; // Avoid singular endpoints

                const double r1_norm = std::sqrt(r1_sq);
                const double r2_norm = std::sqrt(r2_sq);
                const double r1_r2 = r1_norm * r2_norm;
                const double dot_r1_r2 = r1.dot(r2);

                const Vec3 cross_r1_r2 = r1.cross(r2);
                const double cross_sq = cross_r1_r2.norm_squared();
                if (cross_sq < 1e-12) continue; // Collinear

                // K = Gamma/4pi * (|r1| + |r2|) / ( |r1||r2|(|r1||r2| + r1·r2) + smoothing )
                const double numerator = (r1_norm + r2_norm);
                const double denominator = r1_r2 * (r1_r2 + dot_r1_r2) + smoothing_term[f];
                const double K = gamma_over_4pi[f] * (numerator / denominator);

                total_vel_at_p.x += K * cross_r1_r2.x;
                total_vel_at_p.y += K * cross_r1_r2.y;
                total_vel_at_p.z += K * cross_r1_r2.z;
            }
            inducedVelocities[p_idx] = total_vel_at_p;
        }
        double t_calc_end = omp_get_wtime();

        if (simParams.logPerf && workload > 1000000) { // Only log significant workloads
            Logger::log("PERF", "Biot-Savart: Prep=" + std::to_string(t_prep_end - t_prep_start) + 
                        "s, Calc=" + std::to_string(t_calc_end - t_calc_start) + 
                        "s, Workload=" + std::to_string(workload) + ", OMP=" + (enableOpenMP?"On":"Off"));
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
        if (simParams.logVerbose) {
            Logger::log("WAKE", "Initializing wake structure at t=0");
        }

        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &currentBladeWake = wake.getBladeWake(0, b); // 获取 t=0, 叶片 b 的数据结构

            // 1. 添加 t=0 的节点 (叶片上的点 + 初始尾迹点)
            for (int i = 0; i < wake.nTrail; ++i) // nTrail = nShed + 1
            {
                Vec3 boundPos = pos.quarterAt(b, 0, i);                                                            
                currentBladeWake.boundNodeIndices[i] = currentBladeWake.addNode({boundPos, initialFreeStreamVel}); 

                Vec3 trailPos = pos.trailAt(b, 0, i);
                currentBladeWake.trailNodeIndices[i] = currentBladeWake.addNode({trailPos, initialFreeStreamVel}); 
            }

            // 2. 计算 t=0 的涡线强度 (基于初始 BEM 结果)
            std::vector<double> gamma_bound(wake.nShed);
            std::vector<double> gamma_trail(wake.nTrail);

            for (int i = 0; i < wake.nShed; ++i) 
            {
                double cl = perf.clAt(b, 0, i);
                double chord = geom.chordShedding[i]; 
                gamma_bound[i] = 0.5 * turbineParams.windSpeed * chord * cl;
            }

            gamma_trail[0] = gamma_bound[0]; 
            for (int i = 1; i < wake.nShed; ++i)
            {
                gamma_trail[i] = gamma_bound[i] - gamma_bound[i - 1];
            }
            gamma_trail[wake.nShed] = -gamma_bound[wake.nShed - 1]; 

            // 3. 添加 t=0 的涡线
            currentBladeWake.boundLineIndices.resize(wake.nShed, -1);
            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.boundNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i + 1];
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_bound[i], VortexLineType::Bound, gamma_bound[i], false, i});
                currentBladeWake.boundLineIndices[i] = lineIdx;
            }

            currentBladeWake.trailingLineIndices.resize(wake.nTrail, -1);
            for (int i = 0; i < wake.nTrail; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.boundNodeIndices[i];
                int lineIdx = currentBladeWake.addLine({startIdx, endIdx, gamma_trail[i], VortexLineType::Trailing, gamma_trail[i], false, i});
                currentBladeWake.trailingLineIndices[i] = lineIdx;
            }

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
        }

        if (simParams.logVerbose) {
            Logger::log("VEL", "Computing initial velocities for t=0");
        }
        UpdateWakeVelocities(wake, turbineParams, 0, geom, simParams);
    }

    void UpdateWakeVelocities(Wake &wake, const TurbineParams &turbineParams, int timestep, const BladeGeometry &geom, const SimParams &simParams)
    {
        double t_start = omp_get_wtime();
        Vec3 initialFreeStreamVel = Vec3(turbineParams.windSpeed, 0.0, 0.0);

        std::vector<Vec3> allNodePositions;
        std::vector<std::pair<int, int>> nodeGlobalIndices; 
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
            return;
        }

        std::vector<Vec3> inducedVelocities;
        computeInducedVelocity(inducedVelocities, allNodePositions, wake, timestep, geom, simParams);

        if (inducedVelocities.size() != allNodePositions.size())
        {
            std::cerr << "Error: Mismatch between number of target points and calculated induced velocities." << std::endl;
            return; 
        }

        for (size_t i = 0; i < inducedVelocities.size(); ++i)
        {
            int bladeIdx = nodeGlobalIndices[i].first;
            int localNodeIdx = nodeGlobalIndices[i].second;
            wake.getNodes(timestep, bladeIdx)[localNodeIdx].velocity = initialFreeStreamVel + inducedVelocities[i];
        }
        double t_end = omp_get_wtime();
        if (simParams.logPerf) {
            Logger::log("PERF", "UpdateWakeVelocities (Convection): " + std::to_string(t_end - t_start) + " s. Targets: " + std::to_string(allNodePositions.size()));
        }
        if (simParams.logVerbose) {
            Logger::log("VEL", "Updated node velocities");
        }
    }

    void AdvanceWakeStructure(Wake &wake,
                              const TurbineParams &turbineParams, const PositionData &pos,
                              const BladeGeometry &geom, const SimParams &simParams, int currentTimestep)
    {
        const double dt = simParams.dt;
        wake.ensureTimeStepExists(currentTimestep);

        Vec3 initialFreeStreamVel = Vec3(turbineParams.windSpeed, 0.0, 0.0);
        std::vector<std::vector<Vec3>> prevPositionsPerBlade(wake.nBlades);
        std::vector<std::vector<Vec3>> prevVelocitiesPerBlade(wake.nBlades);
        std::vector<int> convectedCountPerBlade(wake.nBlades, 0);

        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);
            const BladeWake &prevBladeWake = wake.getBladeWake(currentTimestep - 1, b);

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

            currentBladeWake.nodes.clear();
            currentBladeWake.lines.clear();
            currentBladeWake.boundNodeIndices.assign(wake.nTrail, -1); 
            currentBladeWake.trailNodeIndices.assign(wake.nTrail, -1); 

            currentBladeWake.boundLineIndices.assign(wake.nShed, -1);
            currentBladeWake.trailingLineIndices.assign(wake.nTrail, -1);
            currentBladeWake.shedLineIndices.assign(wake.nShed, -1);

            const int nConvected = static_cast<int>(prevBladeWake.nodes.size());
            convectedCountPerBlade[b] = nConvected;
            prevPositionsPerBlade[b].reserve(nConvected);
            prevVelocitiesPerBlade[b].reserve(nConvected);

            for (size_t i = 0; i < prevBladeWake.nodes.size(); ++i)
            {
                const VortexNode &node = prevBladeWake.nodes[i];
                prevPositionsPerBlade[b].push_back(node.position);
                prevVelocitiesPerBlade[b].push_back(node.velocity);
                Vec3 newPos = node.position + node.velocity * dt;
                currentBladeWake.addNode({newPos, initialFreeStreamVel});
            }

            for (int i = 0; i < wake.nTrail; ++i)
            {
                currentBladeWake.trailNodeIndices[i] = prevBladeWake.boundNodeIndices[i];

                Vec3 boundPos = pos.quarterAt(b, currentTimestep, i);
                int newIdx = currentBladeWake.addNode({boundPos, initialFreeStreamVel});
                currentBladeWake.boundNodeIndices[i] = newIdx;
            }

            for (const auto &line : prevBladeWake.lines)
            {
                if (line.type == VortexLineType::Trailing || line.type == VortexLineType::Shed)
                {
                    currentBladeWake.addLine({line.startNodeIdx, line.endNodeIdx, line.gamma, line.type, line.initial_gamma, line.in_far_wake, line.segment_index});
                }
            }

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

            for (int i = 0; i < wake.nShed; ++i)
            {
                int startIdx = currentBladeWake.trailNodeIndices[i];
                int endIdx = currentBladeWake.trailNodeIndices[i + 1];
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

        // Explicit predictor-corrector (Heun):
        // x* = x^n + dt u^n
        // x^{n+1} = x^n + 0.5*dt*(u^n + u*)
        if (simParams.timeScheme == TimeMarchingScheme::PredictorCorrector)
        {
            std::vector<Vec3> allNodePositions;
            std::vector<std::pair<int, int>> nodeGlobalIndices;
            for (int b = 0; b < wake.nBlades; ++b)
            {
                const auto &nodes = wake.getNodes(currentTimestep, b);
                for (size_t n_idx = 0; n_idx < nodes.size(); ++n_idx)
                {
                    allNodePositions.push_back(nodes[n_idx].position);
                    nodeGlobalIndices.push_back({b, static_cast<int>(n_idx)});
                }
            }

            if (!allNodePositions.empty())
            {
                std::vector<Vec3> inducedVelocities;
                computeInducedVelocity(inducedVelocities, allNodePositions, wake, currentTimestep, geom, simParams);
                for (size_t i = 0; i < inducedVelocities.size(); ++i)
                {
                    const int b = nodeGlobalIndices[i].first;
                    const int localNodeIdx = nodeGlobalIndices[i].second;
                    const int nConvected = convectedCountPerBlade[b];
                    if (localNodeIdx >= nConvected) continue; // Newly added bound nodes are fixed by kinematics.

                    const Vec3 uStar = initialFreeStreamVel + inducedVelocities[i];
                    const Vec3 &xN = prevPositionsPerBlade[b][localNodeIdx];
                    const Vec3 &uN = prevVelocitiesPerBlade[b][localNodeIdx];
                    wake.getNodes(currentTimestep, b)[localNodeIdx].position = xN + (uN + uStar) * (0.5 * dt);
                }
            }
        }

    }

    // Kutta循环 update vortex strength
    void kuttaJoukowskiIteration(Wake &wake, PerformanceData &perf, const BladeGeometry &geom, NodeAxes &axes,
                                 const PositionData &pos, VelBCS &velBCS,
                                 std::vector<AirfoilData> &airfoils, const TurbineParams &turbineParams, const SimParams &simParams)
    {
        const int maxIterations = simParams.kuttaMaxIterations;
        const double convergenceThreshold = simParams.kuttaTolerance;
        const double relaxationFactor = simParams.kuttaRelaxation;
        bool if_verbose = false;

        int currentTimestep = wake.bladeWakes.size() - 1;
        if (currentTimestep < 0) return;

        std::vector<Vec3> controlPointPositions;
        std::vector<std::pair<int, int>> controlPointMapping;
        controlPointPositions.reserve(wake.nBlades * wake.nShed);

        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nShed; ++i)
            {
                try {
                    controlPointPositions.push_back(pos.boundAt(b, currentTimestep, i));
                    controlPointMapping.push_back({b, i});
                } catch(...) {}
            }
        }
        if (controlPointPositions.empty()) return;

        size_t splitIdx = std::numeric_limits<size_t>::max();
        {
            const BladeWake &bw = wake.getBladeWake(currentTimestep, 0);
            size_t min_idx = bw.lines.size();
            
            auto update_min = [&](const std::vector<int>& indices) {
                for(int idx : indices) {
                    if (idx >= 0 && static_cast<size_t>(idx) < min_idx) min_idx = idx;
                }
            };
            
            update_min(bw.boundLineIndices);
            update_min(bw.trailingLineIndices);
            update_min(bw.shedLineIndices);
            
            splitIdx = min_idx;
            
            if (simParams.logVerbose) {
                Logger::log("DEBUG", "Wake Split Analysis at t=" + std::to_string(currentTimestep));
                Logger::log("DEBUG", "  Total Lines: " + std::to_string(bw.lines.size()));
                Logger::log("DEBUG", "  Split Index: " + std::to_string(splitIdx));
                Logger::log("DEBUG", "  Static Lines (0-" + std::to_string(splitIdx) + "): " + std::to_string(splitIdx));
                Logger::log("DEBUG", "  Dynamic Lines (" + std::to_string(splitIdx) + "-" + std::to_string(bw.lines.size()) + "): " + std::to_string(bw.lines.size() - splitIdx));
                Logger::log("DEBUG", "  Targets: " + std::to_string(controlPointPositions.size()));
            }
        }

        // --- 2. Compute Static Induction (ONCE) ---
        // Influence of everything BEFORE the first dynamic line.
        double t_static_start = omp_get_wtime();
        std::vector<Vec3> v_static(controlPointPositions.size());
        computeInducedVelocity(v_static, controlPointPositions, wake, currentTimestep, geom, simParams, 0, splitIdx);
        double t_static_end = omp_get_wtime();
        if (simParams.logVerbose) {
            Logger::log("DEBUG", "Static Induction Time: " + std::to_string(t_static_end - t_static_start) + " s");
        }

        // --- 3. Iteration Loop ---
        double t_loop_start = omp_get_wtime();
        std::vector<Vec3> v_dynamic(controlPointPositions.size());
        std::vector<std::vector<double>> updatedBoundGammas(wake.nBlades, std::vector<double>(wake.nShed));
        
        double max_dg = 0.0;
        for (int iter = 0; iter < maxIterations; ++iter)
        {
            if (if_verbose) std::cout << "  Kutta Iter " << iter << std::endl;
            max_dg = 0.0;

            computeInducedVelocity(v_dynamic, controlPointPositions, wake, currentTimestep, geom, simParams, splitIdx, std::numeric_limits<size_t>::max());

            for (size_t k = 0; k < controlPointPositions.size(); ++k)
            {
                int b = controlPointMapping[k].first;
                int i = controlPointMapping[k].second;

                Vec3 uind = v_static[k] + v_dynamic[k];
                
                const Vec3 &bxn = axes.bxnAt(b, currentTimestep, i); 
                const Vec3 &byn = axes.bynAt(b, currentTimestep, i); 
                const Vec3 &bzn = axes.bznAt(b, currentTimestep, i); 

                Vec3 uind_blade_frame;
                uind_blade_frame.x = bxn.dot(uind); 
                uind_blade_frame.y = byn.dot(uind); 
                uind_blade_frame.z = bzn.dot(uind); 

                Vec3 vel_eff_blade_frame = velBCS.at(b, currentTimestep, i) + uind_blade_frame;
                perf.setRelativeVelocityAt(b, currentTimestep, i) = vel_eff_blade_frame;

                double Vinf_eff_squared = vel_eff_blade_frame.norm_squared();
                double Vinf_eff = std::sqrt(Vinf_eff_squared);
                double aoa_deg = std::atan2(-vel_eff_blade_frame.y, vel_eff_blade_frame.x) * 180.0 / M_PI;

                int airfoilIdx = geom.airfoilIndex[i];
                double cl_value = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cl, aoa_deg);
                double cd_value = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cd, aoa_deg);

                double chord = geom.chordShedding[i];
                double phi = 0.0;
                double tip_factor = 1.0;
                if (simParams.tipLossModel != TipLossModel::Off &&
                    compute_phi(b, currentTimestep, i, pos, bxn, byn, bzn, vel_eff_blade_frame, phi))
                {
                    tip_factor = std::clamp(
                        compute_tip_loss_factor(geom.rShedding[i], phi, turbineParams, simParams.tipLossModel),
                        0.0, 1.0);
                }
                const double cl_effective = cl_value * tip_factor;

                perf.setAoaAt(b, currentTimestep, i) = aoa_deg;
                perf.setClAt(b, currentTimestep, i) = cl_effective;
                perf.setCdAt(b, currentTimestep, i) = cd_value;

                // Compute rotor-normal/tangential force per unit span (N/m)
                double vxy = std::sqrt(vel_eff_blade_frame.x * vel_eff_blade_frame.x +
                                       vel_eff_blade_frame.y * vel_eff_blade_frame.y);
                if (vxy > 1e-12)
                {
                    double ex = vel_eff_blade_frame.x / vxy;
                    double ey = vel_eff_blade_frame.y / vxy;

                    double L = 0.5 * turbineParams.rho * vxy * vxy * chord * cl_effective;
                    double D = 0.5 * turbineParams.rho * vxy * vxy * chord * cd_value;

                    // Force in blade frame (x: chord, y: normal)
                    double fx_b = L * (-ey) + D * (-ex);
                    double fy_b = L * (ex) + D * (-ey);

                    Vec3 force_ics = bxn * fx_b + byn * fy_b;

                    Vec3 axis_ics(-1.0, 0.0, 0.0); // rotor thrust axis (opposite wind +x)
                    Vec3 radial_ics = pos.boundAt(b, currentTimestep, i) - pos.hubAt(currentTimestep);
                    radial_ics.x = 0.0; // project to rotor plane (y-z)
                    double radial_norm = radial_ics.norm();
                    if (radial_norm > 1e-12) radial_ics = radial_ics * (1.0 / radial_norm);

                    Vec3 tangential_ics = axis_ics.cross(radial_ics);
                    double tan_norm = tangential_ics.norm();
                    if (tan_norm > 1e-12) tangential_ics = tangential_ics * (1.0 / tan_norm);

                    perf.setFnAt(b, currentTimestep, i) = force_ics.dot(axis_ics);
                    perf.setFtAt(b, currentTimestep, i) = force_ics.dot(tangential_ics);
                }
                else
                {
                    perf.setFnAt(b, currentTimestep, i) = 0.0;
                    perf.setFtAt(b, currentTimestep, i) = 0.0;
                }

                double gamma_required = 0.5 * Vinf_eff * chord * cl_effective;

                BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);
                int boundLineIdx = currentBladeWake.boundLineIndices[i];
                VortexLine &boundLine = currentBladeWake.lines[boundLineIdx]; 

                double gamma_current = boundLine.gamma;
                double dg = gamma_required - gamma_current;
                boundLine.gamma = gamma_current + relaxationFactor * dg;

                if (static_cast<size_t>(i) < updatedBoundGammas[b].size()) updatedBoundGammas[b][i] = boundLine.gamma;
                
                max_dg = std::max(max_dg, std::abs(dg));
            }

            for (int b = 0; b < wake.nBlades; ++b)
            {
                BladeWake &currentBladeWake = wake.getBladeWake(currentTimestep, b);

                for (int i = 0; i < wake.nTrail; ++i) 
                {
                    if (static_cast<size_t>(i) >= currentBladeWake.trailingLineIndices.size()) continue;
                    int idx = currentBladeWake.trailingLineIndices[i];
                    if (idx < 0) continue;
                    
                    VortexLine &line = currentBladeWake.lines[idx];
                    double g = 0.0;
                    if (i == 0) {
                        if (updatedBoundGammas[b].size() > 0) g = updatedBoundGammas[b][0];
                    } else if (i < wake.nShed) {
                        g = updatedBoundGammas[b][i] - updatedBoundGammas[b][i-1];
                    } else { // Last one
                        g = -updatedBoundGammas[b][wake.nShed-1];
                    }
                    line.gamma = g;
                }

                for (int i = 0; i < wake.nShed; ++i)
                {
                    if (static_cast<size_t>(i) >= currentBladeWake.shedLineIndices.size()) continue;
                    int idx = currentBladeWake.shedLineIndices[i];
                    if (idx < 0) continue;

                    VortexLine &line = currentBladeWake.lines[idx];
                    double prev_g = (static_cast<size_t>(i) < currentBladeWake.prevGammaBound.size()) ? currentBladeWake.prevGammaBound[i] : 0.0;
                    double curr_g = updatedBoundGammas[b][i];
                    line.gamma = prev_g - curr_g;
                }
            }

            if (max_dg < convergenceThreshold)
            {
                if (simParams.logVerbose) {
                    Logger::log("KUTTA", "Converged: " + std::to_string(iter + 1) + " iters, |dGamma| = " + std::to_string(max_dg));
                }
                break;
            }
        }
        
        double t_loop_end = omp_get_wtime();
        if (simParams.logVerbose) {
            Logger::log("DEBUG", "Dynamic Loop Time: " + std::to_string(t_loop_end - t_loop_start) + " s (" + std::to_string(controlPointPositions.size()) + " per iter)");
        }

        if (max_dg >= convergenceThreshold)
        {
             Logger::log("KUTTA", "WARNING: Failed to converge. Max |dGamma| = " + std::to_string(max_dg));
        }

        for (int b = 0; b < wake.nBlades; ++b)
        {
            for (int i = 0; i < wake.nShed; ++i)
            {
                perf.setBoundGammaAt(b, currentTimestep, i) = updatedBoundGammas[b][i];
            }
        }
        
    }

    void ApplyGammaDecayAndRemoval(Wake &wake, int timestep, const TurbineParams &turbineParams, const SimParams &simParams)
    {
        const double x_far_wake_D = 5.0; 
        const double x_far_wake_m = x_far_wake_D * turbineParams.rTip * 2.0;

        const double alpha = 0.2;                  
        const double gamma_threshold_ratio = 0.01; 

        for (int b = 0; b < wake.nBlades; ++b)
        {
            BladeWake &bladeWake = wake.getBladeWake(timestep, b);

            for (VortexLine &line : bladeWake.lines)
            {
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
                        line.in_far_wake = true; 
                    }
                }
            }

            auto &lines = bladeWake.lines;
            auto original_size = lines.size();

            lines.erase(
                std::remove_if(lines.begin(), lines.end(),
                               [gamma_threshold_ratio](const VortexLine &line)
                               {
                                   if (!line.in_far_wake)
                                       return false;

                                   if (std::abs(line.initial_gamma) < 1e-9)
                                   {
                                       return std::abs(line.gamma) < 1e-7;
                                   }

                                   return std::abs(line.gamma / line.initial_gamma) < gamma_threshold_ratio;
                               }),
                lines.end());

            if (original_size > lines.size())
            {
                // std::cout << "  - Blade " << b << ": Removed " << original_size - lines.size() << " weak vortices." << std::endl;
            }
        }
    }

} // namespace fvw
