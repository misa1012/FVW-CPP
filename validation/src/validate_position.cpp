// Example: ./fvw_cpp --validate --section=position --blade=0 --timestep=0 --nodeType="trailing"

#include "../include/validate_position.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <algorithm>

namespace fvw
{

    void validatePosition(const PositionData &pos, const SimParams &simParams,
                          const TurbineParams &turbineParams, const BladeGeometry &geom,
                          int blade, int timestep, const std::string &nodeType, bool saveToFile)
    {
        // 设置输出精度
        std::cout << std::fixed << std::setprecision(6);

        // 验证时间步有效性
        if (timestep < 0 || timestep >= simParams.timesteps)
        {
            std::cerr << "Error: Invalid timestep (" << timestep << "), must be in [0, "
                      << simParams.timesteps - 1 << "]" << std::endl;
            return;
        }

        // 确定验证的叶片
        int startBlade = (blade == -1) ? 0 : blade;
        int endBlade = (blade == -1) ? turbineParams.nBlades : blade + 1;
        if (blade >= turbineParams.nBlades || blade < -1)
        {
            std::cerr << "Error: Invalid blade (" << blade << "), must be in [-1, "
                      << turbineParams.nBlades - 1 << "]" << std::endl;
            return;
        }

        // 验证节点类型
        std::string nodeTypeLower = nodeType;
        std::transform(nodeTypeLower.begin(), nodeTypeLower.end(), nodeTypeLower.begin(), ::tolower);
        bool validateTrailing = (nodeTypeLower == "all" || nodeTypeLower == "trailing");
        bool validateShedding = (nodeTypeLower == "all" || nodeTypeLower == "shedding");
        if (!validateTrailing && !validateShedding)
        {
            std::cerr << "Error: Invalid nodeType (" << nodeType << "), must be 'all', 'trailing', or 'shedding'" << std::endl;
            return;
        }

        // 1. 验证 Hub
        if (nodeTypeLower == "all")
        {
            std::cout << "\nValidating hub at t=" << timestep << ":" << std::endl;
            Vec3 hub_actual = pos.hubAt(timestep);
            Vec3 hub_expected(0.0, 0.0, 90.0);
            std::cout << "  Expected: (" << hub_expected.x << ", " << hub_expected.y << ", "
                      << hub_expected.z << ")" << std::endl;
            std::cout << "  Actual: (" << hub_actual.x << ", " << hub_actual.y << ", "
                      << hub_actual.z << ")" << std::endl;
            if (std::abs(hub_actual.x - hub_expected.x) > 1e-6 ||
                std::abs(hub_actual.y - hub_expected.y) > 1e-6 ||
                std::abs(hub_actual.z - hub_expected.z) > 1e-6)
            {
                std::cerr << "Error: Hub position incorrect at t=" << timestep << std::endl;
            }

            // 2. 验证 Platform
            std::cout << "\nValidating platform at t=" << timestep << ":" << std::endl;
            Vec3 platform_actual = pos.platformAt(timestep);
            Vec3 platform_expected(0.0, 0.0, 0.0);
            std::cout << "  Expected: (" << platform_expected.x << ", " << platform_expected.y << ", "
                      << platform_expected.z << ")" << std::endl;
            std::cout << "  Actual: (" << platform_actual.x << ", " << platform_actual.y << ", "
                      << platform_actual.z << ")" << std::endl;
            if (std::abs(platform_actual.x - platform_expected.x) > 1e-6 ||
                std::abs(platform_actual.y - platform_expected.y) > 1e-6 ||
                std::abs(platform_actual.z - platform_expected.z) > 1e-6)
            {
                std::cerr << "Error: Platform position incorrect at t=" << timestep << std::endl;
            }
        }

        // 3. 验证 Trailing 节点
        if (validateTrailing)
        {
            for (int b = startBlade; b < endBlade; ++b)
            {
                std::cout << "\nValidating trailing nodes for blade " << b << " at t=" << timestep
                          << " (lead, quarter, trail):" << std::endl;
                if (geom.rTrailing.size() != static_cast<size_t>(turbineParams.nSegments + 1))
                {
                    std::cerr << "Error: rTrailing size (" << geom.rTrailing.size()
                              << ") != nSegments + 1 (" << turbineParams.nSegments + 1 << ")" << std::endl;
                }
                for (size_t i = 0; i < geom.rTrailing.size(); ++i)
                {
                    std::cout << "i=" << i << ":" << std::endl;
                    Vec3 lead = pos.leadAt(b, timestep, i);
                    Vec3 quarter = pos.quarterAt(b, timestep, i);
                    Vec3 trail = pos.trailAt(b, timestep, i);
                    std::cout << "  lead: (" << lead.x << ", " << lead.y << ", " << lead.z << ")" << std::endl;
                    std::cout << "  quarter: (" << quarter.x << ", " << quarter.y << ", " << quarter.z << ")" << std::endl;
                    std::cout << "  trail: (" << trail.x << ", " << trail.y << ", " << trail.z << ")" << std::endl;
                }
            }
        }

        // 4. 验证 Shedding 节点
        if (validateShedding)
        {
            for (int b = startBlade; b < endBlade; ++b)
            {
                std::cout << "\nValidating shedding nodes for blade " << b << " at t=" << timestep
                          << " (colloc, bound, end):" << std::endl;
                if (geom.rShedding.size() != static_cast<size_t>(turbineParams.nSegments))
                {
                    std::cerr << "Error: rShedding size (" << geom.rShedding.size()
                              << ") != nSegments (" << turbineParams.nSegments << ")" << std::endl;
                }
                for (size_t i = 0; i < geom.rShedding.size(); ++i)
                {
                    std::cout << "i=" << i << ":" << std::endl;
                    Vec3 colloc = pos.collocAt(b, timestep, i);
                    Vec3 bound = pos.boundAt(b, timestep, i);
                    Vec3 end = pos.endAt(b, timestep, i);
                    std::cout << "  colloc: (" << colloc.x << ", " << colloc.y << ", " << colloc.z << ")" << std::endl;
                    std::cout << "  bound: (" << bound.x << ", " << bound.y << ", " << bound.z << ")" << std::endl;
                    std::cout << "  end: (" << end.x << ", " << end.y << ", " << end.z << ")" << std::endl;
                }
            }
        }

        // 保存到 CSV
        if (saveToFile)
        {
            std::ofstream out("position.csv");
            out << std::fixed << std::setprecision(6);
            out << "blade,timestep,node_type,i,x,y,z\n";
            // Hub 和 Platform
            if (nodeTypeLower == "all")
            {
                Vec3 hub_actual = pos.hubAt(timestep);
                Vec3 platform_actual = pos.platformAt(timestep);
                out << "-1," << timestep << ",hub,-1," << hub_actual.x << "," << hub_actual.y << ","
                    << hub_actual.z << "\n";
                out << "-1," << timestep << ",platform,-1," << platform_actual.x << "," << platform_actual.y << ","
                    << platform_actual.z << "\n";
            }
            // Trailing 节点
            if (validateTrailing)
            {
                for (int b = startBlade; b < endBlade; ++b)
                {
                    for (size_t i = 0; i < geom.rTrailing.size(); ++i)
                    {
                        Vec3 lead = pos.leadAt(b, timestep, i);
                        Vec3 quarter = pos.quarterAt(b, timestep, i);
                        Vec3 trail = pos.trailAt(b, timestep, i);
                        out << b << "," << timestep << ",lead," << i << "," << lead.x << "," << lead.y << ","
                            << lead.z << "\n";
                        out << b << "," << timestep << ",quarter," << i << "," << quarter.x << "," << quarter.y << ","
                            << quarter.z << "\n";
                        out << b << "," << timestep << ",trail," << i << "," << trail.x << "," << trail.y << ","
                            << trail.z << "\n";
                    }
                }
            }
            // Shedding 节点
            if (validateShedding)
            {
                for (int b = startBlade; b < endBlade; ++b)
                {
                    for (size_t i = 0; i < geom.rShedding.size(); ++i)
                    {
                        Vec3 colloc = pos.collocAt(b, timestep, i);
                        Vec3 bound = pos.boundAt(b, timestep, i);
                        Vec3 end = pos.endAt(b, timestep, i);
                        out << b << "," << timestep << ",colloc," << i << "," << colloc.x << "," << colloc.y << ","
                            << colloc.z << "\n";
                        out << b << "," << timestep << ",bound," << i << "," << bound.x << "," << bound.y << ","
                            << bound.z << "\n";
                        out << b << "," << timestep << ",end," << i << "," << end.x << "," << end.y << ","
                            << end.z << "\n";
                    }
                }
            }
            out.close();
            std::cout << "\nPosition data saved to position.csv" << std::endl;
        }
    }
} // namespace fvw