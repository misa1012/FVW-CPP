#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "wake.h"
#include "postprocess.h"
#include "geometry.h"

// 新增一个函数，用于将计算出的速度场写入VTK文件
void write_field_to_vtk(const std::string &filename,
                        const std::vector<fvw::Vec3> &grid_points,
                        const std::vector<fvw::Vec3> &velocities,
                        int Nx, int Ny, int Nz);

int main()
{
    const std::string h5_filepath = "/home/shug8104/sa/vortex/postprocess/20250612_3000steps/wake.h5";
    const std::string vtk_filepath = "/home/shug8104/sa/vortex/postprocess/20250612_3000steps/velocity_field_t3000.vtk";
    const int timestep_to_process = 3000;
    const double TURBINE_DIAMETER = 126.0;

    // =======================================================================

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

    // 3. 定义欧拉网格 (与之前相同)
    const double cell_size_m = 5.0; // 每个网格单元的边长（米）
    std::cout << "根据 " << cell_size_m << "m 的分辨率创建欧拉网格..." << std::endl;

    const double x_start_D = -1.0, x_end_D = 10.0;
    const double y_range_D = 1.5;
    const double z_range_D = 1.5;

    const double x_range_m = (x_end_D - x_start_D) * TURBINE_DIAMETER;
    const double y_range_m = y_range_D * TURBINE_DIAMETER;
    const double z_range_m = z_range_D * TURBINE_DIAMETER;

    // 使用 std::ceil 来向上取整，确保网格能完整覆盖指定范围
    const int Nx = static_cast<int>(std::ceil(x_range_m / cell_size_m)) + 1;
    const int Ny = static_cast<int>(std::ceil(y_range_m / cell_size_m)) + 1;
    const int Nz = static_cast<int>(std::ceil(z_range_m / cell_size_m)) + 1;

    std::cout << "正在创建欧拉网格点..." << std::endl;
    std::vector<fvw::Vec3> grid_points;
    grid_points.reserve(Nx * Ny * Nz);

    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                double x = (x_start_D + (x_end_D - x_start_D) * i / (Nx - 1)) * TURBINE_DIAMETER;
                double y = (-y_range_D + (2.0 * y_range_D) * j / (Ny - 1)) * (TURBINE_DIAMETER / 2.0);
                double z = (0.0 + (2.0 * z_range_D) * k / (Nz - 1)) * (TURBINE_DIAMETER / 2.0);
                grid_points.push_back({x, y, z});
            }
        }
    }
    std::cout << "网格创建完毕。维度: " << Nx << " x " << Ny << " x " << Nz
              << ", 总点数: " << grid_points.size() << std::endl;

    // 4. 调用函数，计算网格上每个点的诱导速度 (与之前相同)
    std::cout << "开始计算网格点上的诱导速度 ..." << std::endl;
    auto calc_start = std::chrono::high_resolution_clock::now();

    std::vector<fvw::Vec3> grid_velocities;
    double cutOff = 0.005;
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

// --- VTK写入函数 (与之前相同) ---
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