// 在下游定义几个点，计算其随时间变化的induced velocity

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip> // 用于设置输出精度

// 包含项目中已有的头文件
#include "wake.h"
#include "postprocess.h" // 需要我们刚刚添加的新函数
#include "geometry.h"    // TurbineParams在这里定义

int main()
{
    // 自定义参数
    std::string h5_filepath = "/home/shug8104/sa/vortex/postprocess/20250717_01_chord_base_cutoff_study/0_1/wake.h5";
    std::string csv_filepath = "/home/shug8104/sa/vortex/postprocess/20250717_01_chord_base_cutoff_study/0_1/probe_output.csv";
    double D = 126.0;

    fvw::SimParams simParams;
    simParams.coreType = fvw::VortexCoreType::ChordBasedCore;
    simParams.cutoffParam = 0.1; // Biot-Savart的cutOff参数
    simParams.vortexModel = fvw::VortexModelType::Constant;

    // 你的数据是每10步输出一次，共3000步
    int start_step = 0;
    int end_step = 1000;
    int step_interval = 10;
    double dt = 0.06;

    // 这里需要一个 TurbineParams 对象来传递给 read_wake_snapshot
    // 最好的方法是从HDF5的config组里读取，这里为了简化，我们先手动创建

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

    // 2. 定义你的“虚拟探针”位置
    // 对于每个下游X位置，我们将计算4个点并取平均
    std::vector<double> x_locations_D = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    double offset = 50.0; // 偏移量，根据您的需求设定

    // 3. 准备输出文件
    std::ofstream outfile(csv_filepath);
    outfile << "Timestep,Time,DownstreamX_D,Avg_U_ind,Avg_V_ind,Avg_W_ind\n"; // 修改CSV头部，因为现在输出的是平均值
    outfile << std::fixed << std::setprecision(6); // 设置输出精度

    // 4. 循环遍历HDF5中的所有时间步

    // 创建一个Wake对象，它将在循环中被重复填充
    fvw::Wake wake(turbineParams.nBlades, turbineParams.nSegments, turbineParams.nSegments + 1);

    for (int t = start_step; t <= end_step; t += step_interval)
    {
        std::cout << "正在处理时间步: " << t << "..." << std::endl;

        // 5. 从HDF5中读取数据，填充到wake对象在时间步t的状态中
        fvw::read_wake_snapshot(wake, h5_filepath, t, turbineParams);

        // 对于每个下游X位置
        for (double x_D : x_locations_D)
        {
            std::vector<fvw::Vec3> current_probe_points;
            // 添加四个探针点
            current_probe_points.push_back({x_D * D, 0.0, 90.0 + offset});  // y=0, z=90+offset
            current_probe_points.push_back({x_D * D, 0.0, 90.0 - offset});  // y=0, z=90-offset
            current_probe_points.push_back({x_D * D, offset, 90.0});        // y=offset, z=90
            current_probe_points.push_back({x_D * D, -offset, 90.0});       // y=-offset, z=90

            // 计算这四个探针点的诱导速度
            std::vector<fvw::Vec3> induced_velocities_at_current_probes;
            fvw::computeInducedVelocity(induced_velocities_at_current_probes, current_probe_points, wake, t, turbineParams, geom, simParams);

            // 计算平均诱导速度
            fvw::Vec3 sum_vel = {0.0, 0.0, 0.0};
            for (const auto& vel : induced_velocities_at_current_probes)
            {
                sum_vel.x += vel.x;
                sum_vel.y += vel.y;
                sum_vel.z += vel.z;
            }
            fvw::Vec3 avg_vel = {
                sum_vel.x / current_probe_points.size(),
                sum_vel.y / current_probe_points.size(),
                sum_vel.z / current_probe_points.size()
            };

            // 将平均结果写入CSV文件
            outfile << t << "," << t * dt << "," << x_D << ","
                    << avg_vel.x << "," << avg_vel.y << "," << avg_vel.z << "\n";
        }
    }

    outfile.close();
    std::cout << "\n后处理完成，结果已保存至 " << csv_filepath << std::endl;

    return 0;
}
