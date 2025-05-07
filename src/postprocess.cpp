#include "postprocess.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <map>
#include <vector>
#include <iostream>
#include <H5Cpp.h>

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

    void writeWakeToHDF5(const Wake &wake, const TurbineParams &turbineParams,
                         const std::string &outputFile, int timestep)
    {
        try
        {
            // 打开或创建 HDF5 文件（追加模式）
            H5::H5File file(outputFile, timestep == 0 ? H5F_ACC_TRUNC : H5F_ACC_RDWR);

            // 创建时间步组
            std::string groupName = "/timestep_" + std::to_string(timestep);
            H5::Group timestepGroup = file.createGroup(groupName);

            for (int b = 0; b < wake.nBlades; ++b)
            {
                const BladeWake &bladeWake = wake.getBladeWake(timestep, b);
                const auto &nodes = bladeWake.nodes;
                const auto &lines = bladeWake.lines;

                // 创建叶片组
                std::string bladeGroupName = groupName + "/blade_" + std::to_string(b);
                H5::Group bladeGroup = file.createGroup(bladeGroupName);

                // 1. 写入节点数据
                if (!nodes.empty())
                {
                    hsize_t nodeDims[2] = {nodes.size(), 6}; // [nNodes, (x,y,z,u,v,w)]
                    H5::DataSpace nodeSpace(2, nodeDims);

                    H5::DataSet nodeDataset = bladeGroup.createDataSet(
                        "nodes", H5::PredType::NATIVE_DOUBLE, nodeSpace);

                    std::vector<double> nodeData(nodes.size() * 6);
                    for (size_t i = 0; i < nodes.size(); ++i)
                    {
                        nodeData[i * 6 + 0] = nodes[i].position.x;
                        nodeData[i * 6 + 1] = nodes[i].position.y;
                        nodeData[i * 6 + 2] = nodes[i].position.z;
                        nodeData[i * 6 + 3] = nodes[i].velocity.x;
                        nodeData[i * 6 + 4] = nodes[i].velocity.y;
                        nodeData[i * 6 + 5] = nodes[i].velocity.z;
                    }

                    nodeDataset.write(&nodeData[0], H5::PredType::NATIVE_DOUBLE);
                }

                // 2. 写入涡量线数据
                if (!lines.empty())
                {
                    hsize_t lineDims[2] = {lines.size(), 4}; // [nLines, (start_idx,end_idx,gamma,type)]
                    H5::DataSpace lineSpace(2, lineDims);

                    H5::DataSet lineDataset = bladeGroup.createDataSet(
                        "lines", H5::PredType::NATIVE_DOUBLE, lineSpace);

                    std::vector<double> lineData(lines.size() * 4);
                    for (size_t i = 0; i < lines.size(); ++i)
                    {
                        lineData[i * 4 + 0] = static_cast<double>(lines[i].startNodeIdx);
                        lineData[i * 4 + 1] = static_cast<double>(lines[i].endNodeIdx);
                        lineData[i * 4 + 2] = lines[i].gamma;
                        lineData[i * 4 + 3] = static_cast<double>(lines[i].type);
                    }

                    lineDataset.write(&lineData[0], H5::PredType::NATIVE_DOUBLE);
                }
            }

            std::cout << "Wrote HDF5 data for timestep " << timestep << " to " << outputFile << std::endl;
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
        }
    }
} // namespace fvw