// Example: ./fvw_cpp --validate --section=position --blade=0 --timestep=0 --nodeType="trailing"
// --validate: 验证所有（geometry 和 position）。
// --section=geometry: 仅验证 geometry。
// --section=position: 仅验证 position。
// --blade=N: 指定叶片（position）。
// --timestep=T: 指定时间步（position）。
// --nodeType=type: 指定节点类型（position）。
// --save-csv: 保存 CSV 文件。

// 解析验证参数
// 目前已从main.cpp中删去
// auto valParams = fvw::parseValidationArgs(argc, argv);
// 验证 geometry（如果需要）
// if (!valParams.section.empty() &&
//     (valParams.section == "all" || valParams.section == "geometry"))
// {
//     fvw::ValidationParams geomParams = valParams;
//     geomParams.section = "geometry"; // 仅验证 geometry
//     fvw::validate(geomParams, simParams, turbineParams, geom, pos);
// }

// // 验证 position（如果需要）
// if (!valParams.section.empty() &&
//     (valParams.section == "all" || valParams.section == "position"))
// {
//     fvw::ValidationParams posParams = valParams;
//     posParams.section = "position"; // 仅验证 position
//     fvw::validate(posParams, simParams, turbineParams, geom, pos);
// }

#include "validate.h"
#include "validate_geometry.h"
#include "validate_position.h"
#include <algorithm>
#include <iostream>

namespace fvw
{

    // 解析命令行的validation参数
    ValidationParams parseValidationArgs(int argc, char *argv[])
    {
        ValidationParams params;
        params.section = ""; // 默认不验证

        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "--validate")
            {
                params.section = "all"; // 默认验证所有部分
            }
            else if (arg.find("--section=") == 0)
            {
                params.section = arg.substr(10);
            }
            else if (arg.find("--blade=") == 0)
            {
                try
                {
                    params.blade = std::stoi(arg.substr(8));
                }
                catch (...)
                {
                    std::cerr << "Warning: Invalid blade value, using default (-1)" << std::endl;
                }
            }
            else if (arg.find("--timestep=") == 0)
            {
                try
                {
                    params.timestep = std::stoi(arg.substr(11));
                }
                catch (...)
                {
                    std::cerr << "Warning: Invalid timestep value, using default (0)" << std::endl;
                }
            }
            else if (arg.find("--nodeType=") == 0)
            {
                params.nodeType = arg.substr(11);
            }
            else if (arg == "--save-csv")
            {
                params.saveToFile = true;
            }
        }

        return params;
    }

    void validate(const ValidationParams &params, const SimParams &simParams,
                  const TurbineParams &turbineParams, const BladeGeometry &geom,
                  const PositionData &pos)
    {
        if (params.section.empty())
        {
            return; // 无验证请求
        }

        std::string sectionLower = params.section;
        std::transform(sectionLower.begin(), sectionLower.end(), sectionLower.begin(), ::tolower);

        bool if_validateGeometry = (sectionLower == "all" || sectionLower == "geometry");
        bool if_validatePosition = (sectionLower == "all" || sectionLower == "position");

        if (!if_validateGeometry && !if_validatePosition)
        {
            std::cerr << "Error: Invalid section (" << params.section
                      << "), must be 'all', 'geometry', or 'position'" << std::endl;
            return;
        }

        if (if_validateGeometry)
        {
            std::cout << "\n=== Validating Geometry ===" << std::endl;
            validateBladeGeometry(geom, turbineParams, params.saveToFile);
        }

        if (if_validatePosition)
        {
            std::cout << "\n=== Validating Position ===" << std::endl;
            validatePosition(pos, simParams, turbineParams, geom, params.blade,
                             params.timestep, params.nodeType, params.saveToFile);
        }
    }

} // namespace fvw