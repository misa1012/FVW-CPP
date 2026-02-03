#include "validate_geometry.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace fvw
{

    void validateBladeGeometry(const BladeGeometry &geom, const TurbineParams &params, bool saveToFile)
    {
        // 设置输出精度
        std::cout << std::fixed << std::setprecision(6);

        // 1. 验证 trailing 节点
        std::cout << "\nValidating trailing nodes (rTrailing, chordTrailing, twistTrailing):" << std::endl;
        if (geom.rTrailing.size() != static_cast<size_t>(params.nSegments + 1))
        {
            std::cerr << "Error: rTrailing size (" << geom.rTrailing.size() << ") != nSegments + 1 ("
                      << params.nSegments + 1 << ")" << std::endl;
        }
        for (size_t i = 0; i < geom.rTrailing.size(); ++i)
        {
            std::cout << "i=" << i << ": r=" << geom.rTrailing[i]
                      << ", chord=" << geom.chordTrailing[i]
                      << ", twist=" << geom.twistTrailing[i] << " deg" << std::endl;
            // 检查单调性
            if (i > 0 && geom.rTrailing[i] <= geom.rTrailing[i - 1])
            {
                std::cerr << "Error: rTrailing not monotonic at i=" << i << std::endl;
            }
        }
        // 检查边界
        if (std::abs(geom.rTrailing[0] - params.rHub) > 1e-6 ||
            std::abs(geom.rTrailing.back() - params.rTip) > 1e-6)
        {
            std::cerr << "Error: rTrailing bounds incorrect: ["
                      << geom.rTrailing[0] << ", " << geom.rTrailing.back() << "]" << std::endl;
        }

        // 2. 验证 shedding 节点
        std::cout << "\nValidating shedding nodes (rShedding, chordShedding, twistShedding, airfoilIndex):" << std::endl;
        if (geom.rShedding.size() != static_cast<size_t>(params.nSegments))
        {
            std::cerr << "Error: rShedding size (" << geom.rShedding.size() << ") != nSegments ("
                      << params.nSegments << ")" << std::endl;
        }
        for (size_t i = 0; i < geom.rShedding.size(); ++i)
        {
            std::cout << "i=" << i << ": r=" << geom.rShedding[i]
                      << ", chord=" << geom.chordShedding[i]
                      << ", twist=" << geom.twistShedding[i]
                      << ", airfoil=" << geom.airfoilIndex[i] << std::endl;
            // 检查 chord 范围
            if (geom.chordShedding[i] <= 0)
            {
                std::cerr << "Error: chordShedding non-positive at i=" << i << std::endl;
            }
        }

        // 3. 验证节点位置
        std::cout << "\nValidating node positions (colloc, bound, end):" << std::endl;
        for (size_t i = 0; i < geom.colloc.size(); ++i)
        {
            std::cout << "i=" << i << ":\n";
            std::cout << "  colloc: (" << geom.colloc[i].x << ", " << geom.colloc[i].y << ", " << geom.colloc[i].z << ")\n";
            std::cout << "  bound: (" << geom.bound[i].x << ", " << geom.bound[i].y << ", " << geom.bound[i].z << ")\n";
            std::cout << "  end: (" << geom.end[i].x << ", " << geom.end[i].y << ", " << geom.end[i].z << ")\n";
            // 检查 colloc.x = 0.25 * chord
            double expectedCollocX = 0.25 * geom.chordShedding[i];
            if (std::abs(geom.colloc[i].x - expectedCollocX) > 1e-6)
            {
                std::cerr << "Error: colloc.x incorrect at i=" << i << std::endl;
            }
            // 检查 z = rShedding
            if (std::abs(geom.colloc[i].z - geom.rShedding[i]) > 1e-6 ||
                std::abs(geom.bound[i].z - geom.rShedding[i]) > 1e-6 ||
                std::abs(geom.end[i].z - geom.rShedding[i]) > 1e-6)
            {
                std::cerr << "Error: z-coordinate mismatch at i=" << i << std::endl;
            }
        }

        // Trailing 节点 (lead, quarter, trail)
        std::cout << "\nValidating trailing node positions (lead, quarter, trail):" << std::endl;
        for (size_t i = 0; i < geom.trail.size(); ++i)
        {
            std::cout << "Trailing i=" << i << ":\n";
            std::cout << "  lead: (" << geom.lead[i].x << ", " << geom.lead[i].y << ", " << geom.lead[i].z << ")\n";
            std::cout << "  quarter: (" << geom.quarter[i].x << ", " << geom.quarter[i].y << ", " << geom.quarter[i].z << ")\n";
            std::cout << "  trail: (" << geom.trail[i].x << ", " << geom.trail[i].y << ", " << geom.trail[i].z << ")\n";
            // 检查 z = rTrailing
            if (std::abs(geom.lead[i].z - geom.rTrailing[i]) > 1e-6 ||
                std::abs(geom.quarter[i].z - geom.rTrailing[i]) > 1e-6 ||
                std::abs(geom.trail[i].z - geom.rTrailing[i]) > 1e-6)
            {
                std::cerr << "Error: z-coordinate mismatch for trailing at i=" << i << std::endl;
            }
        }

        // 4. 保存到 CSV
        if (saveToFile)
        {
            std::ofstream out("geometry.csv");
            out << std::fixed << std::setprecision(6);
            out << "node_type,i,r,chord,twist,airfoil_index,x,y,z\n";
            // Trailing nodes
            for (size_t i = 0; i < geom.rTrailing.size(); ++i)
            {
                out << "trailing," << i << "," << geom.rTrailing[i] << "," << geom.chordTrailing[i]
                    << "," << geom.twistTrailing[i] << ",-1,"
                    << geom.trail[i].x << "," << geom.trail[i].y << "," << geom.trail[i].z << "\n";
                out << "trailing_lead," << i << "," << geom.rTrailing[i] << "," << geom.chordTrailing[i]
                    << "," << geom.twistTrailing[i] << ",-1,"
                    << geom.lead[i].x << "," << geom.lead[i].y << "," << geom.lead[i].z << "\n";
                out << "trailing_quarter," << i << "," << geom.rTrailing[i] << "," << geom.chordTrailing[i]
                    << "," << geom.twistTrailing[i] << ",-1,"
                    << geom.quarter[i].x << "," << geom.quarter[i].y << "," << geom.quarter[i].z << "\n";
            }
            // Shedding nodes
            for (size_t i = 0; i < geom.rShedding.size(); ++i)
            {
                out << "shedding," << i << "," << geom.rShedding[i] << "," << geom.chordShedding[i]
                    << "," << geom.twistShedding[i] << "," << geom.airfoilIndex[i]
                    << "," << geom.colloc[i].x << "," << geom.colloc[i].y << "," << geom.colloc[i].z << "\n";
                out << "shedding_bound," << i << "," << geom.rShedding[i] << "," << geom.chordShedding[i]
                    << "," << geom.twistShedding[i] << "," << geom.airfoilIndex[i]
                    << "," << geom.bound[i].x << "," << geom.bound[i].y << "," << geom.bound[i].z << "\n";
                out << "shedding_end," << i << "," << geom.rShedding[i] << "," << geom.chordShedding[i]
                    << "," << geom.twistShedding[i] << "," << geom.airfoilIndex[i]
                    << "," << geom.end[i].x << "," << geom.end[i].y << "," << geom.end[i].z << "\n";
            }
            out.close();
            std::cout << "\nGeometry data saved to geometry.csv" << std::endl;
        }
    }

} // namespace fvw