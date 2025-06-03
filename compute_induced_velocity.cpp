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
#include <chrono>

void saveInducedVelocityToHDF5(const std::vector<std::vector<fvw::Vec3>>& vel_induced_results,
                               int num_points, int timesteps, const std::string& filename) {
    try {
        H5::H5File file(filename, H5F_ACC_TRUNC);
        hsize_t dims[2] = {static_cast<hsize_t>(timesteps), 3};
        H5::DataSpace dataspace(2, dims);

        for (int p = 0; p < num_points; ++p) {
            std::string dataset_name = "induced_velocity_point_" + std::to_string(p);
            H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> data(timesteps * 3);
            for (int t = 0; t < timesteps; ++t) {
                data[t * 3 + 0] = vel_induced_results[t][p].x;
                data[t * 3 + 1] = vel_induced_results[t][p].y;
                data[t * 3 + 2] = vel_induced_results[t][p].z;
            }
            dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
        }
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
    }
}

int main()
{
    // 总计时开始
    auto total_start = std::chrono::high_resolution_clock::now();

    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

    // Section 1. Define parameters
    // Call geometry.h
    // Simulation parameters
    fvw::SimParams simParams;
    simParams.dt = 0.06;
    simParams.totalTime = 3.0;
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
    auto airfoils = fvw::readAirfoils("/home/sa/vortex/FVW-CPP/data/NREL_files/");
    std::cout << "Airfoil profiles read-in completed." << std::endl;

    // Section 3. Calculate the initial position
    // Call position.h
    // Compute positions
    fvw::PositionData pos(turbineParams.nBlades, simParams.timesteps,
                          turbineParams.nSegments + 1, turbineParams.nSegments);
    fvw::computePositions(pos, simParams, turbineParams, geom);
    std::cout << "Position computation completed." << std::endl;

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
        fvw::writeWakeToHDF5(wake, pos, perf, velICS, velBCS, turbineParams, "../results/wake.h5", 0);
        fvw::writeConfigToHDF5(geom, turbineParams, simParams, "../results/wake.h5");
    }

    // 对特定的点计算induced velocity
    // 定义观察点
    const double D = 126.0;
    const double observe_z = 90.0 + 50.0;
    std::vector<fvw::Vec3> observe_points = {
        fvw::Vec3(0.5 * D, 0.0, observe_z),
        fvw::Vec3(1.0 * D, 0.0, observe_z),
        fvw::Vec3(1.5 * D, 0.0, observe_z),
        fvw::Vec3(2.0 * D, 0.0, observe_z)};
    const int num_points = observe_points.size();
    // 计算 results 数据集
    std::vector<std::vector<fvw::Vec3>> vel_induced_results(simParams.timesteps, std::vector<fvw::Vec3>(num_points));

    // --- 主时步推进循环 ---
    for (int t = 1; t < simParams.timesteps; ++t)
    {
        auto step_start = std::chrono::high_resolution_clock::now();

        std::cout << "\n--- Advancing timestep " << t << "/" << simParams.timesteps - 1
                  << " (physical time: " << t * simParams.dt << " s) ---" << std::endl;

        // 1. 推进尾迹结构到时间步 t
        // 这将对流 t-1 的节点并添加 t 的新附着节点。
        // 注意：AdvanceWakeStructure 已经确保时间步 t 存在。
        fvw::AdvanceWakeStructure(wake, geom, turbineParams, pos, simParams.dt, t);

        // 2. 对时间步 t 执行 Kutta-Joukowski 迭代
        // 这将更新时间步 t 的附着涡线、脱落涡线和分离涡线的 gamma 值。
        fvw::kuttaJoukowskiIteration(wake, perf, geom, axes, turbineParams, pos, velBCS, airfoils);

        // 3. 更新时间步 t 的尾迹节点速度
        // 根据更新后的 gamma 值计算所有节点的总速度（诱导速度 + 自由来流速度）
        fvw::UpdateWakeVelocities(wake, turbineParams, t);

        // 4. 写入 VTK 文件
        fvw::writeWakeToVTK(wake, turbineParams, "../results/output", t);
        fvw::writeWakeToHDF5(wake, pos, perf, velICS, velBCS, turbineParams, "../results/wake.h5", t);

        // 5. 对某些点计算induced velocity并输出hdf5文件
        fvw::computeInducedVelocity(vel_induced_results[t], observe_points, wake, t, turbineParams, 0.001);
        
        auto step_end = std::chrono::high_resolution_clock::now();
        auto step_duration = std::chrono::duration_cast<std::chrono::microseconds>(step_end - step_start);
        std::cout << "[Timing] Timestep " << t << ": "
                  << step_duration.count() / 1e6 << " s" << std::endl;
    }

    // 保存诱导速度到 HDF5（单一文件，多个数据集）
    saveInducedVelocityToHDF5(vel_induced_results, num_points, simParams.timesteps, "../results/induced_velocity.h5");
    std::cout << "Saved induced velocities to ../results/induced_velocity.h5" << std::endl;

    std::cout << "\n[END] Wake computation completed." << std::endl;

    // 总计时结束
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(total_end - total_start);
    std::cout << "[Timing] Total time: "
              << total_duration.count() / 1e6 << " s" << std::endl;

    return 0;
}