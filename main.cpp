#include "airfoil.h"
#include "geometry.h"
#include "position.h"
#include "velocity.h"
#include "performance.h"
#include "bem.h"
#include "wake.h"
#include "postprocess.h"
#include "validate.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <cmath>

void write_field_to_vtk(const std::string &filename,
                        const std::vector<fvw::Vec3> &grid_points,
                        const std::vector<fvw::Vec3> &velocities,
                        int Nx, int Ny, int Nz);

void projectWakeToEulerianGrid(const fvw::Wake &wake,
                               const fvw::TurbineParams &turbineParams,
                               const fvw::SimParams &simParams,
                               const fvw::BladeGeometry &geom,
                               int final_timestep,
                               const std::string &outputPath);

int main(int argc, char *argv[])
{
#ifdef NDEBUG
    std::cout << "This is a Release build." << std::endl;
#else
    std::cout << "This is a Debug build." << std::endl;
#endif

    // --- 控制开关与路径设置 ---
    const bool projectToGrid = false;                                                                               // 是否在模拟结束后将最终尾流投影到欧拉网格
    const std::string outputPath = "/home/shug8104/sa/vortex/postprocess/20250727_01_chord_base_cutoff_study/0_1"; // 设置输出文件目录

    std::filesystem::create_directories(outputPath); // 确保输出目录存在

    // 总计时开始
    auto total_start = std::chrono::high_resolution_clock::now();

    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

    // Section 1. Define parameters
    // Call geometry.h
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

    // Simulation parameters
    fvw::SimParams simParams;
    simParams.dt = 0.06;
    simParams.totalTime = 180.0;
    simParams.timesteps = static_cast<int>(simParams.totalTime / simParams.dt) + 1;
    simParams.outputFrequency = 10; // 每10步输出一次HDF5

    simParams.coreType = fvw::VortexCoreType::ChordBasedCore;
    simParams.cutoffParam = 0.1;                            // 用于控制cutoff的参数，van Garrel就是delta，chordbase就是选弦长的比例
    simParams.vortexModel = fvw::VortexModelType::Constant; // 选择vortex diffusion model

    // --- 扰动实验 ---
    simParams.perturbation.type = fvw::PerturbationType::None;
    simParams.perturbation.amplitude_deg = 0.1; // 扰动幅值：度
    // 将扰动频率设置为与涡轮旋转频率相关的值
    double rotational_freq_hz = turbineParams.omega / (2.0 * M_PI);
    simParams.perturbation.frequency_hz = 4.5 * rotational_freq_hz; // 例如，扰动频率 = 4.5*涡轮旋转频率

    // Compute blade geometry
    auto geom = fvw::computeBladeGeometry(turbineParams);
    std::cout << "Blade geometry computed." << std::endl;

    // Section 2. Read in airfoil files
    // Call airfoil.h
    // Read airfoils
    auto airfoils = fvw::readAirfoils("/home/shug8104/sa/vortex/FVW-CPP/data/NREL_files/");
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
    fvw::InitializeWakeStructure(wake, geom, perf, turbineParams, pos, simParams);

    // timestep=0：写入配置和初始数据

    // 不写入vtk文件以节省时间
    // fvw::writeWakeToVTK(wake, turbineParams, "../results/output", 0);
    std::string h5_filepath = (std::filesystem::path(outputPath) / "wake.h5").string();
    fvw::writeWakeToHDF5(wake, pos, perf, velICS, velBCS, turbineParams, h5_filepath, 0);
    fvw::writeConfigToHDF5(geom, turbineParams, simParams, h5_filepath);

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
        fvw::kuttaJoukowskiIteration(wake, perf, geom, axes, turbineParams, pos, velBCS, airfoils, simParams);

        // *** 在这里调用新的衰减和移除函数 ***
        if (simParams.vortexModel == fvw::VortexModelType::GammaDecay)
        {
            std::cout << "Applying gamma decay model for timestep " << t << "..." << std::endl;
            fvw::ApplyGammaDecayAndRemoval(wake, t, turbineParams, simParams);
        }

        // 3. 更新时间步 t 的尾迹节点速度
        // 根据更新后的 gamma 值计算所有节点的总速度（诱导速度 + 自由来流速度）
        fvw::UpdateWakeVelocities(wake, turbineParams, t, geom, simParams);

        // 4. 写入 VTK 文件
        // fvw::writeWakeToVTK(wake, turbineParams, "../results/output", t);
        // fvw::writeWakeToHDF5(wake, pos, perf, velICS, velBCS, turbineParams, "../results/wake.h5", t);

        // --- 按频率输出HDF5 ---
        if (t % simParams.outputFrequency == 0 || t == simParams.timesteps - 1)
        {
            fvw::writeWakeToHDF5(wake, pos, perf, velICS, velBCS, turbineParams, h5_filepath, t);
        }

        auto step_end = std::chrono::high_resolution_clock::now();
        auto step_duration = std::chrono::duration_cast<std::chrono::microseconds>(step_end - step_start);
        std::cout << "[Timing] Timestep " << t << ": "
                  << step_duration.count() / 1e6 << " s" << std::endl;
    }

    // 最后一步输出vtk文件
    // fvw::writeWakeToVTK(wake, turbineParams, "../results/", simParams.timesteps - 1);
    std::cout << "\n[END] Wake computation completed." << std::endl;

    // 总计时结束
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(total_end - total_start);
    std::cout << "[Timing] Total time: "
              << total_duration.count() / 1e6 << " s" << std::endl;

    // --- 将最终尾流投影到欧拉网格 ---
    if (projectToGrid)
    {
        auto projection_start = std::chrono::high_resolution_clock::now();
        std::cout << "\n--- Projecting final wake to Eulerian grid ---" << std::endl;

        projectWakeToEulerianGrid(wake, turbineParams, simParams, geom, simParams.timesteps - 1, outputPath);

        auto projection_end = std::chrono::high_resolution_clock::now();
        auto projection_duration = std::chrono::duration_cast<std::chrono::microseconds>(projection_end - projection_start);
        std::cout << "[Timing] Projection to Eulerian grid: "
                  << projection_duration.count() / 1e6 << " s" << std::endl;
    }

    return 0;
}

// --- 函数实现 ---
/**
 * @brief 将最终的拉格朗日尾流数据投影到欧拉网格上，并保存为VTK文件。
 * @param wake 最终的尾流数据结构。
 * @param turbineParams 涡轮机参数。
 * @param simParams 模拟参数。
 * @param geom 叶片几何参数。
 * @param final_timestep 要处理的最终时间步索引。
 */
void projectWakeToEulerianGrid(const fvw::Wake &wake,
                               const fvw::TurbineParams &turbineParams,
                               const fvw::SimParams &simParams,
                               const fvw::BladeGeometry &geom,
                               int final_timestep,
                               const std::string &outputPath)
{
    // ===================== 参数设置 =====================
    std::filesystem::path path(outputPath);
    path /= "velocity_field_final_step.vtk";
    const std::string vtk_filepath = path.string();

    const double TURBINE_DIAMETER = turbineParams.rTip * 2.0;
    const double hub_height = 90.0; // 假设轮毂高度，可根据需要修改

    // 网格参数选择
    const bool use_uniform_grid = false; // true: 均匀网格, false: 非均匀网格
    const double res_high_m = 1.0;       // 高分辨率区域的网格大小（米）

    // 网格范围设置
    const double x_start_m = -1.0 * TURBINE_DIAMETER;
    const double x_end_m = 7.0 * TURBINE_DIAMETER;
    const double y_max_m = 1.5 * TURBINE_DIAMETER / 2.0;
    const double z_max_m = 1.5 * TURBINE_DIAMETER / 2.0;

    // 非均匀网格参数 (如果 use_uniform_grid = false)
    const double res_low_m = 5.0;
    const double x_fine_start_m = 0.0 * TURBINE_DIAMETER;
    const double x_fine_end_m = 3.0 * TURBINE_DIAMETER;
    const double yz_fine_range_m = 1.5 * TURBINE_DIAMETER / 2.0;

    // ===================== 网格构建 =====================
    std::cout << "Creating Eulerian grid (" << (use_uniform_grid ? "uniform" : "non-uniform") << ")..." << std::endl;
    std::vector<double> x_coords, y_coords, z_coords;

    if (use_uniform_grid)
    {
        for (double x = x_start_m; x <= x_end_m; x += res_high_m)
            x_coords.push_back(x);
        for (double y = -y_max_m; y <= y_max_m; y += res_high_m)
            y_coords.push_back(y);
        for (double z = hub_height - z_max_m; z <= hub_height + z_max_m; z += res_high_m)
            z_coords.push_back(z);
    }
    else
    {
        // X coordinates
        for (double x = x_start_m; x < x_fine_start_m; x += res_low_m)
            x_coords.push_back(x);
        for (double x = x_fine_start_m; x < x_fine_end_m; x += res_high_m)
            x_coords.push_back(x);
        for (double x = x_fine_end_m; x <= x_end_m; x += res_low_m)
            x_coords.push_back(x);
        // Y coordinates
        for (double y = -y_max_m; y < -yz_fine_range_m; y += res_low_m)
            y_coords.push_back(y);
        for (double y = -yz_fine_range_m; y <= yz_fine_range_m; y += res_high_m)
            y_coords.push_back(y);
        for (double y = yz_fine_range_m + res_low_m; y <= y_max_m; y += res_low_m)
            y_coords.push_back(y);
        // Z coordinates
        for (double z = hub_height - z_max_m; z < hub_height - yz_fine_range_m; z += res_low_m)
            z_coords.push_back(z);
        for (double z = hub_height - yz_fine_range_m; z <= hub_height + yz_fine_range_m; z += res_high_m)
            z_coords.push_back(z);
        for (double z = hub_height + yz_fine_range_m + res_low_m; z <= hub_height + z_max_m; z += res_low_m)
            z_coords.push_back(z);
    }

    const int Nx = x_coords.size();
    const int Ny = y_coords.size();
    const int Nz = z_coords.size();

    std::vector<fvw::Vec3> grid_points;
    grid_points.reserve(Nx * Ny * Nz);
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                grid_points.push_back({x_coords[i], y_coords[j], z_coords[k]});

    std::cout << "Grid created. Dimensions: " << Nx << " x " << Ny << " x " << Nz
              << ", Total points: " << grid_points.size() << std::endl;

    // ===================== 计算速度 =====================
    std::cout << "Calculating induced velocity on grid points..." << std::endl;
    std::vector<fvw::Vec3> grid_velocities;
    fvw::computeInducedVelocity(grid_velocities, grid_points, wake, final_timestep, turbineParams, geom, simParams);

    // ===================== 输出VTK =====================
    std::cout << "Writing velocity field to VTK file: " << vtk_filepath << std::endl;
    write_field_to_vtk(vtk_filepath, grid_points, grid_velocities, Nx, Ny, Nz);
    std::cout << "Post-processing successfully completed!" << std::endl;
}

/**
 * @brief 将计算出的速度场写入VTK文件，格式为STRUCTURED_GRID。
 * @param filename 输出的VTK文件名。
 * @param grid_points 网格点的坐标集合。
 * @param velocities 每个网格点对应的速度向量集合。
 * @param Nx 网格在X方向的维度。
 * @param Ny 网格在Y方向的维度。
 * @param Nz 网格在Z方向的维度。
 */
void write_field_to_vtk(const std::string &filename,
                        const std::vector<fvw::Vec3> &grid_points,
                        const std::vector<fvw::Vec3> &velocities,
                        int Nx, int Ny, int Nz)
{
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open())
    {
        std::cerr << "Error: Unable to open VTK file " << filename << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Wake Velocity Field\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    vtkFile << "POINTS " << grid_points.size() << " float\n";
    for (const auto &p : grid_points)
        vtkFile << p.x << " " << p.y << " " << p.z << "\n";

    vtkFile << "POINT_DATA " << velocities.size() << "\n";
    vtkFile << "VECTORS Velocity float\n";
    for (const auto &v : velocities)
        vtkFile << v.x << " " << v.y << " " << v.z << "\n";

    vtkFile.close();
}