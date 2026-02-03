#ifndef SIMULATION_RUNNER_H
#define SIMULATION_RUNNER_H

#include <string>
#include <filesystem>
#include <chrono>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <memory>

// FVW Includes
#include "io/config.h"
#include "core/geometry.h"
#include "core/position.h"
#include "core/velocity.h"
#include "models/performance.h"
#include "models/bem.h"
#include "core/wake.h"
#include "io/postprocess.h"
#include "models/airfoil.h"

namespace fvw {

struct RunMetadata {
    std::string config_path;
    std::string config_text;
    std::string exe_path;
    std::string cli_args;
};

class SimulationRunner {
public:
    SimulationRunner(const PerturbationConfig& pc, const GlobalConfig& config, const std::string& rootOutput, const RunMetadata& runMeta);
    
    // 初始化模拟环境，加载资源
    void initialize();
    
    // 运行时间步循环
    void run();
    
    // 执行结束后的清理或额外后处理
    void finalize();

private:
    PerturbationConfig m_pc;
    GlobalConfig m_globalConfig;
    std::string m_rootOutput;
    RunMetadata m_runMeta;
    std::string m_caseOutDir;
    std::string m_h5Filepath;
    std::string m_dataRoot;
    double m_totalSimSeconds = 0.0;

    // Simulation Objects
    TurbineParams m_turbineParams;
    SimParams m_simParams;
    
    // Geometry & Flow State
    BladeDefinition m_bladeDef;
    BladeGeometry m_geom;
    std::string m_geometryPath;
    
    // Airfoils
    // Airfoils
    std::vector<AirfoilData> m_airfoils;

    // State Containers (initialized in initialize())
    // Note: Using unique_ptr or keeping them as members if they are copyable/movable.
    // Given the current design where they are large objects, members are fine but constructor needs care.
    // Or we can use std::optional or unique_ptr for delayed initialization.
    // Let's use unique_ptr for explicit control over lifetime and initialization timing.
    std::unique_ptr<PositionData> m_pos;
    std::unique_ptr<VelICS> m_velICS;
    std::unique_ptr<VelBCS> m_velBCS;
    std::unique_ptr<NodeAxes> m_axes;
    std::unique_ptr<PerformanceData> m_perf;
    std::unique_ptr<Wake> m_wake;

    // Helper functions
    void resolve_paths();
    void load_resources();
    void apply_perturbation();
    std::string reset_case_output(const std::string &root, const std::string &case_name, bool clean = true);
    void write_run_manifest();
};

} // namespace fvw

#endif // SIMULATION_RUNNER_H
