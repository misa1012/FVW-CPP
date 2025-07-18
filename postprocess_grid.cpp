// 用于把某一时刻的三维数据map到欧拉网格上
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cmath>

#include "wake.h"
#include "postprocess.h"
#include "geometry.h"

// VTK写入函数的声明
void write_field_to_vtk(const std::string &filename,
                        const std::vector<fvw::Vec3> &grid_points,
                        const std::vector<fvw::Vec3> &velocities,
                        int Nx, int Ny, int Nz);

int main()
{
    // ===================== 参数设置 =====================
    const std::string h5_filepath = "/home/shug8104/sa/vortex/postprocess/20250727_01_chord_base_cutoff_study/0_1/wake.h5";
    const std::string vtk_filepath = "/home/shug8104/sa/vortex/postprocess/20250727_01_chord_base_cutoff_study/0_1/velocity_field_chord_t1000_0_1.vtk";
    const int timestep_to_process = 1000;
    const double TURBINE_DIAMETER = 126.0;
    const double hub_height = 90.0;

    // 涡尾流和涡机参数
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

    // ===================== 数据读取 =====================
    fvw::Wake wake_snapshot(turbineParams.nBlades, turbineParams.nSegments, turbineParams.nSegments + 1);
    std::cout << "正在从 " << h5_filepath << " 读取时间步 " << timestep_to_process << " 的尾流状态..." << std::endl;
    fvw::read_wake_snapshot(wake_snapshot, h5_filepath, timestep_to_process, turbineParams);

    if (wake_snapshot.bladeWakes.empty() || wake_snapshot.getBladeWake(timestep_to_process, 0).nodes.empty()) {
        std::cerr << "错误: 未能加载或尾流为空，程序终止。" << std::endl;
        return 1;
    }
    std::cout << "尾流状态重建成功。" << std::endl;

    // ===================== 网格构建 =====================
    std::cout << "正在创建欧拉网格点 (" << (use_uniform_grid ? "均匀" : "非均匀") << ")..." << std::endl;
    std::vector<double> x_coords, y_coords, z_coords;

    if (use_uniform_grid) {
        for (double x = x_start_m; x <= x_end_m; x += res_high_m) x_coords.push_back(x);
        for (double y = -y_max_m; y <= y_max_m; y += res_high_m) y_coords.push_back(y);
        for (double z = hub_height - z_max_m; z <= hub_height + z_max_m; z += res_high_m) z_coords.push_back(z);
    } else {
        for (double x = x_start_m; x < x_fine_start_m; x += res_low_m) x_coords.push_back(x);
        for (double x = x_fine_start_m; x < x_fine_end_m; x += res_high_m) x_coords.push_back(x);
        for (double x = x_fine_end_m; x <= x_end_m; x += res_low_m) x_coords.push_back(x);

        for (double y = -y_max_m; y < -yz_fine_range_m; y += res_low_m) y_coords.push_back(y);
        for (double y = -yz_fine_range_m; y <= yz_fine_range_m; y += res_high_m) y_coords.push_back(y);
        for (double y = yz_fine_range_m + res_low_m; y <= y_max_m; y += res_low_m) y_coords.push_back(y);

        for (double z = 0; z < hub_height - yz_fine_range_m; z += res_low_m) z_coords.push_back(z);
        for (double z = hub_height - yz_fine_range_m; z <= hub_height + yz_fine_range_m; z += res_high_m) z_coords.push_back(z);
        for (double z = hub_height + yz_fine_range_m + res_low_m; z <= z_max_m; z += res_low_m) z_coords.push_back(z);
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

    std::cout << "网格创建完毕。维度: " << Nx << " x " << Ny << " x " << Nz
              << ", 总点数: " << grid_points.size() << std::endl;

    // ===================== 计算速度 =====================
    std::cout << "开始计算网格点上的诱导速度 ..." << std::endl;
    auto calc_start = std::chrono::high_resolution_clock::now();

    std::vector<fvw::Vec3> grid_velocities;
    fvw::SimParams simParams;
    simParams.coreType = fvw::VortexCoreType::VanGarrel;
    simParams.cutoffParam = 0.1;    // Biot-Savart的cutOff参数
    simParams.vortexModel = fvw::VortexModelType::Constant;
    fvw::computeInducedVelocity(grid_velocities, grid_points, wake_snapshot, timestep_to_process, turbineParams, geom, simParams);

    auto calc_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(calc_end - calc_start);
    std::cout << "速度场计算完成，耗时: " << duration.count() << " s." << std::endl;

    // ===================== 输出VTK =====================
    std::cout << "正在将速度场写入VTK文件: " << vtk_filepath << std::endl;
    write_field_to_vtk(vtk_filepath, grid_points, grid_velocities, Nx, Ny, Nz);
    std::cout << "\n后处理成功完成！" << std::endl;

    return 0;
}

// --- VTK写入函数  ---
void write_field_to_vtk(const std::string &filename,
                        const std::vector<fvw::Vec3> &grid_points,
                        const std::vector<fvw::Vec3> &velocities,
                        int Nx, int Ny, int Nz)
{
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        std::cerr << "错误: 无法打开VTK文件 " << filename << std::endl;
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
