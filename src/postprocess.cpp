#include "postprocess.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <map>
#include <vector>
#include <iostream>

namespace fvw
{
    void writeWakeToVTK(const Wake &wake, const TurbineParams &turbineParams,
                        const std::string &outputDir, int timestep)
    {
        // 确保输出目录存在
        std::filesystem::create_directories(outputDir);

        // 构造输出文件名
        std::stringstream filename;
        filename << outputDir << "/wake_timestep_" << timestep << ".vtk";
        std::ofstream vtkFile(filename.str());

        if (!vtkFile.is_open())
        {
            std::cerr << "Error: Cannot open VTK file for writing: " << filename.str() << std::endl;
            return;
        }

        // --- Write VTK Header ---
        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "Wake data for timestep " << timestep << "\n";
        vtkFile << "ASCII\n";
        vtkFile << "DATASET POLYDATA\n";

        // 1. 收集所有节点
        std::vector<Vec3> uniqueNodes;
        std::map<Vec3, int> nodeIndexMap;                 // 位置到索引的映射
        std::vector<std::pair<int, int>> lineConnections; // 每条线的 (startIdx, endIdx)
        std::vector<double> gammaValues;                  // 每条线的 gamma

        for (int b = 0; b < wake.nBlades; ++b)
        {
            const BladeWake &bladeWake = wake.getBladeWake(timestep, b);
            const auto &nodes = bladeWake.nodes;
            const auto &lines = bladeWake.lines;

            // 遍历节点
            for (const auto &node : nodes)
            {
                Vec3 pos = node.position;
                if (nodeIndexMap.find(pos) == nodeIndexMap.end())
                {
                    int newIndex = uniqueNodes.size();
                    uniqueNodes.push_back(pos);
                    nodeIndexMap[pos] = newIndex;
                }
            }

            // 遍历涡量线
            for (const auto &line : lines)
            {
                if (line.startNodeIdx >= nodes.size() || line.endNodeIdx >= nodes.size())
                {
                    std::cerr << "Warning: Invalid node index in line for blade " << b
                              << ", timestep " << timestep << std::endl;
                    continue;
                }
                Vec3 startPos = nodes[line.startNodeIdx].position;
                Vec3 endPos = nodes[line.endNodeIdx].position;
                int startIdx = nodeIndexMap[startPos];
                int endIdx = nodeIndexMap[endPos];
                lineConnections.push_back({startIdx, endIdx});
                gammaValues.push_back(line.gamma);
            }
        }

        // 2. 写入节点坐标
        vtkFile << "POINTS " << uniqueNodes.size() << " float\n";
        for (const auto &node : uniqueNodes)
        {
            vtkFile << node.x << " " << node.y << " " << node.z << "\n";
        }

        // 3. 写入涡量线连接关系
        vtkFile << "LINES " << lineConnections.size() << " " << lineConnections.size() * 3 << "\n";
        for (const auto &conn : lineConnections)
        {
            vtkFile << "2 " << conn.first << " " << conn.second << "\n";
        }

        // 4. 写入涡量强度 (gamma)
        vtkFile << "CELL_DATA " << lineConnections.size() << "\n";
        vtkFile << "SCALARS gamma float 1\n";
        vtkFile << "LOOKUP_TABLE default\n";
        for (double gamma : gammaValues)
        {
            vtkFile << gamma << "\n";
        }

        vtkFile.close();
        std::cout << "Wrote VTK file: " << filename.str() << std::endl;
    }
} // namespace fvw