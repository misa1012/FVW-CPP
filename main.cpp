#include "airfoil.h"
#include "geometry.h"
#include "position.h"
#include "velocity.h"
#include "performance.h"
#include "bem.h"
#include "wake.h"
#include "postprocess.h"
#include <iostream>
#include <filesystem>

#include "validate.h"

int main(int argc, char *argv[])
{
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

    // 解析验证参数
    auto valParams = fvw::parseValidationArgs(argc, argv);

    // Section 1. Define parameters
    // Call geometry.h
    // Simulation parameters
    fvw::SimParams simParams;
    simParams.dt = 0.06;
    simParams.totalTime = 1.2;
    simParams.timesteps = static_cast<int>(simParams.totalTime / simParams.dt) + 1;
    // Turbine parameters
    fvw::TurbineParams turbineParams;
    turbineParams.windSpeed = 11.4;
    turbineParams.rho = 1.23;
    turbineParams.rHub = 1.5;
    turbineParams.rTip = 63.0;
    turbineParams.nBlades = 3;
    turbineParams.nSegments = 18;
    turbineParams.tsr = 7.0;
    turbineParams.omega = turbineParams.tsr * turbineParams.windSpeed / turbineParams.rTip;

    // Compute blade geometry
    auto geom = fvw::computeBladeGeometry(turbineParams);
    std::cout << "Blade geometry computed." << std::endl;

    // Section 2. Read in airfoil files
    // Call airfoil.h
    // Read airfoils
    auto airfoils = fvw::readAirfoils("../data/NREL_files/");
    std::cout << "Airfoil profiles read-in completed." << std::endl;

    // Section 3. Calculate the initial position
    // Call position.h
    // Compute positions
    fvw::PositionData pos(turbineParams.nBlades, simParams.timesteps,
                          turbineParams.nSegments + 1, turbineParams.nSegments);
    fvw::computePositions(pos, simParams, turbineParams, geom);
    std::cout << "Position computation completed." << std::endl;

    // 验证 geometry（如果需要）
    if (!valParams.section.empty() &&
        (valParams.section == "all" || valParams.section == "geometry"))
    {
        fvw::ValidationParams geomParams = valParams;
        geomParams.section = "geometry"; // 仅验证 geometry
        fvw::validate(geomParams, simParams, turbineParams, geom, pos);
    }

    // 验证 position（如果需要）
    if (!valParams.section.empty() &&
        (valParams.section == "all" || valParams.section == "position"))
    {
        fvw::ValidationParams posParams = valParams;
        posParams.section = "position"; // 仅验证 position
        fvw::validate(posParams, simParams, turbineParams, geom, pos);
    }

    // Call velocity.h
    // Compute velocities
    // Section 4. Compute velocities
    fvw::VelICS velICS(turbineParams.nBlades, simParams.timesteps, turbineParams.nSegments);
    fvw::computeVelICS(velICS, pos, simParams, turbineParams);

    fvw::VelBCS velBCS(turbineParams.nBlades, simParams.timesteps, turbineParams.nSegments);
    fvw::NodeAxes axes(turbineParams.nBlades, simParams.timesteps,
                       turbineParams.nSegments + 1, turbineParams.nSegments);
    fvw::computeVelBCS(velBCS, velICS, axes, pos, simParams, turbineParams);
    std::cout << "Velocity computation completed." << std::endl;

    // Initialize performance data
    fvw::PerformanceData perf(turbineParams.nBlades, simParams.timesteps, turbineParams.nSegments);
    std::cout << "Performance data is initialized." << std::endl;

    // Compute BEM
    fvw::computeBEM(perf, geom, turbineParams, airfoils);
    std::cout << "BEM computation completed." << std::endl;

    // Initialize the wake for the first timestep
    fvw::Wake wake(turbineParams.nBlades, turbineParams.nSegments, turbineParams.nSegments + 1);
    // t=0
    fvw::InitializeWakeStructure(wake, geom, perf, turbineParams, pos, simParams.dt);

    // timestep=0：写入配置和初始数据
    {
        fvw::writeWakeToVTK(wake, turbineParams, "../results/output", 0);
        fvw::writeWakeToHDF5(wake, perf, turbineParams, "../results/wake.h5", 0);
        fvw::writeConfigToHDF5(geom, turbineParams, simParams, "../results/wake.h5");
    }

    // --- 主时步推进循环 ---
    for (int t = 1; t < simParams.timesteps; ++t)
    {
        std::cout << "\n--- Advancing timestep " << t << "/" << simParams.timesteps - 1
                  << " (physical time: " << t * simParams.dt << " s) ---" << std::endl;

        // 1. 推进尾迹结构到时间步 t
        // 这将对流 t-1 的节点并添加 t 的新附着节点。
        // 注意：AdvanceWakeStructure 已经确保时间步 t 存在。
        fvw::AdvanceWakeStructure(wake, geom, perf, turbineParams, pos, simParams.dt, t);

        // 2. 对时间步 t 执行 Kutta-Joukowski 迭代
        // 这将更新时间步 t 的附着涡线、脱落涡线和分离涡线的 gamma 值。
        fvw::kuttaJoukowskiIteration(wake, perf, geom, axes, turbineParams, pos, velBCS, airfoils);

        // 3. 更新时间步 t 的尾迹节点速度
        // 根据更新后的 gamma 值计算所有节点的总速度（诱导速度 + 自由来流速度）。
        fvw::UpdateWakeVelocities(wake, turbineParams, t);

        // 4. 写入 VTK 文件
        fvw::writeWakeToVTK(wake, turbineParams, "../results/output", t);
        fvw::writeWakeToHDF5(wake, perf, turbineParams, "../results/wake.h5", t);
    }

    std::cout << "\n[END] Wake computation completed." << std::endl;

    return 0;
}