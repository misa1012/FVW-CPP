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

    // --- 参数定义 ---
    const std::string h5_filepath = "/home/shug8104/sa/vortex/postprocess/20250612_3000steps/wake.h5";
    const std::string vtk_filepath = "/home/shug8104/sa/vortex/postprocess/20250612_3000steps/velocity_field_t3000_0_1.vtk";
    const int timestep_to_process = 3000;
    const double TURBINE_DIAMETER = 126.0;

    fvw::TurbineParams turbineParams;
    turbineParams.nBlades = 3;
    turbineParams.nSegments = 18;

    // 1. 创建一个空的Wake对象
    fvw::Wake wake_snapshot(turbineParams.nBlades, turbineParams.nSegments, turbineParams.nSegments + 1);

    // 2. 从HDF5中读取数据，填充到我们刚刚创建的wake_snapshot对象中
    std::cout << "正在从 " << h5_filepath << " 读取时间步 " << timestep_to_process << " 的尾流状态..." << std::endl;
    fvw::read_wake_snapshot(wake_snapshot, h5_filepath, timestep_to_process, turbineParams);

    // 检查是否加载成功
    if (wake_snapshot.bladeWakes.empty() || wake_snapshot.getBladeWake(timestep_to_process, 0).nodes.empty())
    {
        std::cerr << "错误: 未能加载或尾流为空，程序终止。" << std::endl;
        return 1;
    }
    std::cout << "尾流状态重建成功。" << std::endl;

    // =======================================================================
    // --- 3. 定义并创建非均匀欧拉网格---
    // =======================================================================
    std::cout << "正在创建非均匀欧拉网格点..." << std::endl;

    // a. 定义各区域的分辨率和边界
    const double res_high_m = 1.0; // 高分辨率区域的单元尺寸 (米)
    const double res_low_m = 5.0;  // 低分辨率区域的单元尺寸 (米)

    // 定义高分辨率“盒子”的范围
    const double x_fine_start_m = 0.0 * TURBINE_DIAMETER;
    const double x_fine_end_m = 3.0 * TURBINE_DIAMETER;

    // 定义整个计算域的范围
    const double x_start_m = -1.0 * TURBINE_DIAMETER;
    const double x_end_m = 3.0 * TURBINE_DIAMETER;
    const double y_max_m = 1.5 * TURBINE_DIAMETER / 2.0; // +/- 1.5R
    const double z_max_m = 1.5 * TURBINE_DIAMETER;       // 0 到 3R

    // b. 分段生成各方向的坐标点
    std::vector<double> x_coords, y_coords, z_coords;

    // 生成 X 坐标 (非均匀)
    for (double x = x_start_m; x < x_fine_start_m; x += res_low_m)
        x_coords.push_back(x);
    for (double x = x_fine_start_m; x < x_fine_end_m; x += res_high_m)
        x_coords.push_back(x);
    for (double x = x_fine_end_m; x <= x_end_m; x += res_low_m)
        x_coords.push_back(x);

    // Y和Z方向也使用同样逻辑，但在中心区域加密
    const double yz_fine_range_m = 1.5 * TURBINE_DIAMETER / 2.0; // 在 +/- 1R 范围内加密
    for (double y = -y_max_m; y < -yz_fine_range_m; y += res_low_m)
        y_coords.push_back(y);
    for (double y = -yz_fine_range_m; y <= yz_fine_range_m; y += res_high_m)
        y_coords.push_back(y);
    for (double y = yz_fine_range_m + res_low_m; y <= y_max_m; y += res_low_m)
        y_coords.push_back(y);

    // Z方向从轮毂中心(假设为90m)附近加密
    const double hub_height = 90.0;
    for (double z = 0; z < hub_height - yz_fine_range_m; z += res_low_m)
        z_coords.push_back(z);
    for (double z = hub_height - yz_fine_range_m; z <= hub_height + yz_fine_range_m; z += res_high_m)
        z_coords.push_back(z);
    for (double z = hub_height + yz_fine_range_m + res_low_m; z <= z_max_m; z += res_low_m)
        z_coords.push_back(z);

    const int Nx = x_coords.size();
    const int Ny = y_coords.size();
    const int Nz = z_coords.size();

    // c. 组合坐标点生成最终的网格点云
    std::vector<fvw::Vec3> grid_points;
    grid_points.reserve(Nx * Ny * Nz);
    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                grid_points.push_back({x_coords[i], y_coords[j], z_coords[k]});
            }
        }
    }

    std::cout << "非均匀网格创建完毕。维度: " << Nx << " x " << Ny << " x " << Nz
              << ", 总点数: " << grid_points.size() << std::endl;
    // =======================================================================

    // --- 4. & 5. 计算速度并写入VTK  ---
    std::cout << "开始计算网格点上的诱导速度 ..." << std::endl;
    auto calc_start = std::chrono::high_resolution_clock::now();

    std::vector<fvw::Vec3> grid_velocities;
    double cutOff = 0.1;
    fvw::computeInducedVelocity(grid_velocities, grid_points, wake_snapshot, timestep_to_process, turbineParams, cutOff);

    auto calc_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(calc_end - calc_start);
    std::cout << "速度场计算完成，耗时: " << duration.count() << " s." << std::endl;

    // 5. 将结果写入VTK文件 (与之前相同)
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
    if (!vtkFile.is_open())
    {
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
    {
        vtkFile << p.x << " " << p.y << " " << p.z << "\n";
    }

    vtkFile << "POINT_DATA " << velocities.size() << "\n";
    vtkFile << "VECTORS Velocity float\n";
    for (const auto &v : velocities)
    {
        vtkFile << v.x << " " << v.y << " " << v.z << "\n";
    }

    vtkFile.close();
}
