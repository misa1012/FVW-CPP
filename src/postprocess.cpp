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

    void writeWakeToHDF5(const Wake &wake, const PerformanceData &perf,
                         const TurbineParams &turbineParams, const std::string &outputFile,
                         int timestep)
    {
        try
        {
            // 打开或创建 HDF5 文件（追加模式）
            H5::H5File file(outputFile, timestep == 0 ? H5F_ACC_TRUNC : H5F_ACC_RDWR);

            // 创建时间步父组
            H5::Group timestepsGroup;
            if (timestep == 0)
            {
                timestepsGroup = file.createGroup("/timesteps");
            }
            else
            {
                timestepsGroup = file.openGroup("/timesteps");
            }

            // 创建时间步组
            std::string groupName = "/timesteps/timestep_" + std::to_string(timestep);
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

                // 3. 写入性能数据 (cl, cd, aoa)
                hsize_t perfDims[2] = {static_cast<hsize_t>(turbineParams.nSegments), 3}; // [nSegments, (cl,cd,aoa)]
                H5::DataSpace perfSpace(2, perfDims);

                H5::DataSet perfDataset = bladeGroup.createDataSet(
                    "perf", H5::PredType::NATIVE_DOUBLE, perfSpace);

                std::vector<double> perfData(turbineParams.nSegments * 3);
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    perfData[i * 3 + 0] = perf.clAt(b, timestep, i);
                    perfData[i * 3 + 1] = perf.cdAt(b, timestep, i);
                    perfData[i * 3 + 2] = perf.aoaAt(b, timestep, i);
                }
                perfDataset.write(&perfData[0], H5::PredType::NATIVE_DOUBLE);
            }

            std::cout << "Wrote HDF5 data for timestep " << timestep << " to " << outputFile << std::endl;
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
        }
    }

    void writeConfigToHDF5(const BladeGeometry &geom, const TurbineParams &turbineParams,
                           const SimParams &simParams, const std::string &outputFile)
    {
        try
        {
             // 打开或创建 HDF5 文件（覆盖模式）
             H5::H5File file(outputFile, H5F_ACC_RDWR);

            // 1. 写入几何数据
            H5::Group geomGroup = file.createGroup("/geometry");

            // 写入 rShedding [nSegments]
            hsize_t rDims[1] = {geom.rShedding.size()};
            H5::DataSpace rSpace(1, rDims);
            H5::DataSet rDataset = geomGroup.createDataSet("r_shed", H5::PredType::NATIVE_DOUBLE, rSpace);
            rDataset.write(geom.rShedding.data(), H5::PredType::NATIVE_DOUBLE);

            // 写入 chordShedding [nSegments]
            hsize_t chordDims[1] = {geom.chordShedding.size()};
            H5::DataSpace chordSpace(1, chordDims);
            H5::DataSet chordDataset = geomGroup.createDataSet("chord", H5::PredType::NATIVE_DOUBLE, chordSpace);
            chordDataset.write(geom.chordShedding.data(), H5::PredType::NATIVE_DOUBLE);

            // 写入 twistShedding [nSegments]
            hsize_t twistDims[1] = {geom.twistShedding.size()};
            H5::DataSpace twistSpace(1, twistDims);
            H5::DataSet twistDataset = geomGroup.createDataSet("twist", H5::PredType::NATIVE_DOUBLE, twistSpace);
            twistDataset.write(geom.twistShedding.data(), H5::PredType::NATIVE_DOUBLE);

            // 写入 airfoil_index [nSegments]
            hsize_t airfoilDims[1] = {geom.airfoilIndex.size()};
            H5::DataSpace airfoilSpace(1, airfoilDims);
            H5::DataSet airfoilDataset = geomGroup.createDataSet("airfoil_index", H5::PredType::NATIVE_INT, airfoilSpace);
            airfoilDataset.write(geom.airfoilIndex.data(), H5::PredType::NATIVE_INT);

            // 2. 写入仿真参数
            H5::Group simGroup = file.createGroup("/simulation");

            // 写入标量参数
            hsize_t scalarDims[1] = {1};
            H5::DataSpace scalarSpace(1, scalarDims);

            H5::DataSet dtDataset = simGroup.createDataSet("dt", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            dtDataset.write(&simParams.dt, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet totalTimeDataset = simGroup.createDataSet("total_time", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            totalTimeDataset.write(&simParams.totalTime, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet nBladesDataset = simGroup.createDataSet("n_blades", H5::PredType::NATIVE_INT, scalarSpace);
            int nBlades = turbineParams.nBlades;
            nBladesDataset.write(&nBlades, H5::PredType::NATIVE_INT);

            H5::DataSet windSpeedDataset = simGroup.createDataSet("wind_speed", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            windSpeedDataset.write(&turbineParams.windSpeed, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet rhoDataset = simGroup.createDataSet("rho", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            rhoDataset.write(&turbineParams.rho, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet rHubDataset = simGroup.createDataSet("r_hub", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            rHubDataset.write(&turbineParams.rHub, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet rTipDataset = simGroup.createDataSet("r_tip", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            rTipDataset.write(&turbineParams.rTip, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet tsrDataset = simGroup.createDataSet("tsr", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            tsrDataset.write(&turbineParams.tsr, H5::PredType::NATIVE_DOUBLE);

            H5::DataSet omegaDataset = simGroup.createDataSet("omega", H5::PredType::NATIVE_DOUBLE, scalarSpace);
            omegaDataset.write(&turbineParams.omega, H5::PredType::NATIVE_DOUBLE);

            // 写入时间步数组
            hsize_t timestepsDims[1] = {static_cast<hsize_t>(simParams.timesteps)};
            H5::DataSpace timestepsSpace(1, timestepsDims);
            H5::DataSet timestepsDataset = simGroup.createDataSet("timesteps", H5::PredType::NATIVE_INT, timestepsSpace);
            std::vector<int> timestepsData(simParams.timesteps);
            for (int i = 0; i < simParams.timesteps; ++i)
            {
                timestepsData[i] = i;
            }
            timestepsDataset.write(timestepsData.data(), H5::PredType::NATIVE_INT);

            std::cout << "Wrote configuration (geometry and simulation parameters) to HDF5 file" << std::endl;
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "HDF5 error in writeConfigToHDF5: " << e.getDetailMsg() << std::endl;
        }
    }
} // namespace fvw