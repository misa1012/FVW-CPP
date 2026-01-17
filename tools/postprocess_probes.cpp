#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip> // 用于设置输出精度
#include <cstdlib> // 用于 atof
#include <filesystem> // C++17 for path manipulation

// 包含项目中已有的头文件
#include "core/wake.h"
#include "io/postprocess.h" // 需要我们刚刚添加的新函数
#include "core/geometry.h"    // TurbineParams在这里定义
#include "io/config.h"

// 引入 std::filesystem 命名空间
namespace fs = std::filesystem;

int main(int argc, char* argv[]) // 修改 main 函数签名以接收命令行参数
{
    // --- 命令行参数处理 ---
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <path_to_wake.h5> <offset_value> [path_to_config.json]" << std::endl;
        std::cerr << "Example: " << argv[0] << " .../wake.h5 50.0" << std::endl;
        return 1; // 退出程序，表示错误
    }

    std::string h5_filepath = argv[1];
    double offset = std::atof(argv[2]); // 将字符串 offset 转换为 double
    
    // Load Config
    std::string config_path = "config.json";
    if (argc >= 4) {
        config_path = argv[3];
    } else {
        if (!fs::exists(config_path) && fs::exists("../config.json")) config_path = "../config.json";
    }

    fvw::GlobalConfig config;
    try {
        config = fvw::ConfigLoader::load(config_path);
    } catch(const std::exception& e) {
        std::cerr << "Error loading config: " << e.what() << "\nUsing defaults.\n";
        return 1;
    }

    // 构造 CSV 输出路径
    fs::path wake_h5_path(h5_filepath);
    std::string output_dir = wake_h5_path.parent_path().string(); // 获取 wake.h5 所在的目录
    std::string csv_filename = "probe_output_" + std::to_string(static_cast<int>(offset)) + ".csv"; // 将 offset 转换为 int 用于文件名，如果你需要浮点数精度，可能需要其他方式格式化
    std::string csv_filepath = output_dir + "/" + csv_filename;

    std::cout << "Wake H5 Path: " << h5_filepath << std::endl;
    std::cout << "Offset Value: " << offset << std::endl;
    std::cout << "CSV Output Path: " << csv_filepath << std::endl;

    // --- 自定义参数 (From Config) ---
    double D = 2.0 * config.turbine.rTip;
    
    // 你的数据是每10步输出一次，共3000步
    // Config should match original logic.
    int start_step = 0;
    int end_step = config.sim.timesteps - 1;
    int step_interval = config.sim.outputFrequency;
    double dt = config.sim.dt;

    // Load geometry (logic similar to others)
    std::string turbine_model = "NREL_5MW"; // should be from config
    std::string data_root = "data/" + turbine_model;
    std::string geometry_path = data_root + "/blade_geometry.csv";
    if (!fs::exists(geometry_path)) {
        data_root = "../data/" + turbine_model;
        geometry_path = data_root + "/blade_geometry.csv";
    }

    // Compute blade geometry
    // Load geometry from file
    auto bladeDef = fvw::loadBladeDefinition(geometry_path);
    if (bladeDef.r.empty()) {
          std::cerr << "Fatal: Could not load blade geometry from " << geometry_path << "\n";
          return 1;
    }
    auto geom = fvw::computeBladeGeometry(config.turbine, bladeDef);
    std::cout << "Blade geometry computed." << std::endl;

    // 2. 定义你的“虚拟探针”位置
    // 对于每个下游X位置，我们将计算4个点并取平均
    std::vector<double> x_locations_D = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    // double offset = 50.0; // 此行现在被命令行参数取代

    // 3. 准备输出文件
    std::ofstream outfile(csv_filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output CSV file: " << csv_filepath << std::endl;
        return 1; // 退出程序，表示错误
    }
    outfile << "Timestep,Time,DownstreamX_D,Avg_U_ind,Avg_V_ind,Avg_W_ind\n"; // 修改CSV头部，因为现在输出的是平均值
    outfile << std::fixed << std::setprecision(6); // 设置输出精度

    // 4. 循环遍历HDF5中的所有时间步
    // 创建一个Wake对象，它将在循环中被重复填充
    fvw::Wake wake(config.turbine.nBlades, config.turbine.nSegments, config.turbine.nSegments + 1);

    for (int t = start_step; t <= end_step; t += step_interval)
    {
        std::cout << "正在处理时间步: " << t << "..." << std::endl;

        // 5. 从HDF5中读取数据，填充到wake对象在时间步t的状态中
        fvw::read_wake_snapshot(wake, h5_filepath, t, config.turbine);

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
            fvw::computeInducedVelocity(induced_velocities_at_current_probes, current_probe_points, wake, t, config.turbine, geom, config.sim);

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