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

    void writeWakeToHDF5(const Wake &wake, const PositionData &pos, const PerformanceData &perf, VelICS &velICS, VelBCS &velBCS,
                         const TurbineParams &turbineParams, const std::string &outputFile,
                         int timestep)
    {
        try
        {
            // 打开或创建 HDF5 文件（追加模式）
            H5::H5File file(outputFile, timestep == 0 ? H5F_ACC_TRUNC : H5F_ACC_RDWR);

            // 1. 创建 Wake 组
            H5::Group wakeGroup;
            if (timestep == 0)
            {
                wakeGroup = file.createGroup("/wake");
            }
            else
            {
                wakeGroup = file.openGroup("/wake");
            }

            std::string wakeTimestepName = "/wake/timestep_" + std::to_string(timestep);
            H5::Group wakeTimestepGroup = wakeGroup.createGroup(wakeTimestepName);

            // 2. 创建 LiftingLine 组
            H5::Group liftingLineGroup;
            if (timestep == 0)
            {
                liftingLineGroup = file.createGroup("/liftingline");
            }
            else
            {
                liftingLineGroup = file.openGroup("/liftingline");
            }

            std::string liftingTimestepName = "/liftingline/timestep_" + std::to_string(timestep);
            H5::Group liftingTimestepGroup = liftingLineGroup.createGroup(liftingTimestepName);

            for (int b = 0; b < wake.nBlades; ++b)
            {
                const BladeWake &bladeWake = wake.getBladeWake(timestep, b);
                const auto &nodes = bladeWake.nodes;
                const auto &lines = bladeWake.lines;

                // Wake 叶片组
                std::string wakeBladeGroupName = wakeTimestepName + "/blade_" + std::to_string(b);
                H5::Group wakeBladeGroup = wakeTimestepGroup.createGroup(wakeBladeGroupName);

                // 1. 写入节点数据
                if (!nodes.empty())
                {
                    hsize_t nodeDims[2] = {nodes.size(), 6}; // [nNodes, (x,y,z,u,v,w)]
                    H5::DataSpace nodeSpace(2, nodeDims);

                    H5::DataSet nodeDataset = wakeBladeGroup.createDataSet("nodes", H5::PredType::NATIVE_DOUBLE, nodeSpace);

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

                    H5::DataSet lineDataset = wakeBladeGroup.createDataSet("lines", H5::PredType::NATIVE_DOUBLE, lineSpace);

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

                // LiftingLine 叶片组
                std::string liftingBladeGroupName = liftingTimestepName + "/blade_" + std::to_string(b);
                H5::Group liftingBladeGroup = liftingTimestepGroup.createGroup(liftingBladeGroupName);

                // 3. 写入性能数据 (cl, cd, aoa)
                hsize_t perfDims[2] = {static_cast<hsize_t>(turbineParams.nSegments), 3}; // [nSegments, (cl,cd,aoa)]
                H5::DataSpace perfSpace(2, perfDims);

                H5::DataSet perfDataset = liftingBladeGroup.createDataSet("perf", H5::PredType::NATIVE_DOUBLE, perfSpace);

                std::vector<double> perfData(turbineParams.nSegments * 3);
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    perfData[i * 3 + 0] = perf.clAt(b, timestep, i);
                    perfData[i * 3 + 1] = perf.cdAt(b, timestep, i);
                    perfData[i * 3 + 2] = perf.aoaAt(b, timestep, i);
                }
                perfDataset.write(&perfData[0], H5::PredType::NATIVE_DOUBLE);

                // 写入relative velocity
                hsize_t velDims[2] = {static_cast<hsize_t>(turbineParams.nSegments), 3};
                H5::DataSpace velSpace(2, velDims);
                H5::DataSet relativeBcsDataset = liftingBladeGroup.createDataSet("relative_velocity_bcs", H5::PredType::NATIVE_DOUBLE, velSpace);

                std::vector<double> relativeBcsData(turbineParams.nSegments * 3);
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    const Vec3 &velocity = perf.relativeVelocityAt(b, timestep, i);
                    relativeBcsData[i * 3 + 0] = velocity.x;
                    relativeBcsData[i * 3 + 1] = velocity.y;
                    relativeBcsData[i * 3 + 2] = velocity.z;
                }
                relativeBcsDataset.write(relativeBcsData.data(), H5::PredType::NATIVE_DOUBLE);

                // 4. 写入诱导速度数据 (u, v, w)
                // H5::Group inducedVelGroup = liftingBladeGroup.createGroup("induced_velocity");

                // hsize_t inducedVelDims[2] = {static_cast<hsize_t>(turbineParams.nSegments), 3}; // [nSegments, (u,v,w)]
                // H5::DataSpace inducedVelSpace(2, inducedVelDims);

                // H5::DataSet bcsDataset = inducedVelGroup.createDataSet("BCS", H5::PredType::NATIVE_DOUBLE, inducedVelSpace);
                // bcsDataset.createAttribute("coordinate_system", H5::StrType(H5::PredType::C_S1, 32), H5::DataSpace(H5S_SCALAR))
                //     .write(H5::StrType(H5::PredType::C_S1, 32), "blade_frame");

                // std::vector<double> bcsData(turbineParams.nSegments * 3);
                // for (int i = 0; i < turbineParams.nSegments; ++i)
                // {
                //     const Vec3 &velocity = perf.inducedVelocityAt(b, timestep, i);
                //     bcsData[i * 3 + 0] = velocity.x; // u
                //     bcsData[i * 3 + 1] = velocity.y; // v
                //     bcsData[i * 3 + 2] = velocity.z; // w
                // }
                // bcsDataset.write(&bcsData[0], H5::PredType::NATIVE_DOUBLE);

                // // ICS (惯性坐标系)
                // H5::DataSet icsDataset = inducedVelGroup.createDataSet("ICS", H5::PredType::NATIVE_DOUBLE, inducedVelSpace);
                // icsDataset.createAttribute("coordinate_system", H5::StrType(H5::PredType::C_S1, 32), H5::DataSpace(H5S_SCALAR))
                //     .write(H5::StrType(H5::PredType::C_S1, 32), "inertial_frame");

                // std::vector<double> icsData(turbineParams.nSegments * 3);
                // for (int i = 0; i < turbineParams.nSegments; ++i)
                // {
                //     const Vec3 &velocity = perf.inducedVelocityICSAt(b, timestep, i);
                //     icsData[i * 3 + 0] = velocity.x; // x
                //     icsData[i * 3 + 1] = velocity.y; // y
                //     icsData[i * 3 + 2] = velocity.z; // z
                // }
                // icsDataset.write(&icsData[0], H5::PredType::NATIVE_DOUBLE);

                // // 2.3 写入 velBCS 数据
                // H5::Group bladeVelocityGroup = liftingBladeGroup.createGroup("blade_velocity");
                // H5::DataSet velBCSDataset = bladeVelocityGroup.createDataSet("velBCS", H5::PredType::NATIVE_DOUBLE, inducedVelSpace);
                // velBCSDataset.createAttribute("coordinate_system", H5::StrType(H5::PredType::C_S1, 32), H5::DataSpace(H5S_SCALAR))
                //     .write(H5::StrType(H5::PredType::C_S1, 32), "blade_frame");

                // std::vector<double> velBCSData(turbineParams.nSegments * 3);
                // for (int i = 0; i < turbineParams.nSegments; ++i)
                // {

                //     const Vec3 &velocity = velBCS.at(b, timestep, i);
                //     velBCSData[i * 3 + 0] = velocity.x; // x
                //     velBCSData[i * 3 + 1] = velocity.y; // y
                //     velBCSData[i * 3 + 2] = velocity.z; // z
                // }
                // velBCSDataset.write(&velBCSData[0], H5::PredType::NATIVE_DOUBLE);

                // H5::DataSet velICSDataset = bladeVelocityGroup.createDataSet("velICS", H5::PredType::NATIVE_DOUBLE, inducedVelSpace);
                // velICSDataset.createAttribute("coordinate_system", H5::StrType(H5::PredType::C_S1, 32), H5::DataSpace(H5S_SCALAR))
                //     .write(H5::StrType(H5::PredType::C_S1, 32), "inertial_frame");

                // std::vector<double> velICSData(turbineParams.nSegments * 3);
                // for (int i = 0; i < turbineParams.nSegments; ++i)
                // {

                //     const Vec3 &velocity = velICS.at(b, timestep, i);
                //     velICSData[i * 3 + 0] = velocity.x; // x
                //     velICSData[i * 3 + 1] = velocity.y; // y
                //     velICSData[i * 3 + 2] = velocity.z; // z
                // }
                // velICSDataset.write(&velICSData[0], H5::PredType::NATIVE_DOUBLE);

                // 2.4 写入 Bound Gamma 数据
                hsize_t gammaDims[2] = {static_cast<hsize_t>(turbineParams.nSegments), 1}; // [nSegments, 1]
                H5::DataSpace gammaSpace(2, gammaDims);
                H5::DataSet gammaDataset = liftingBladeGroup.createDataSet("bound_gamma", H5::PredType::NATIVE_DOUBLE, gammaSpace);

                std::vector<double> gammaData(turbineParams.nSegments);
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    gammaData[i] = perf.boundGammaAt(b, timestep, i);
                }
                gammaDataset.write(&gammaData[0], H5::PredType::NATIVE_DOUBLE);

                // // 2.5 写入 pos.boundAt 数据
                // hsize_t posDims[2] = {static_cast<hsize_t>(turbineParams.nSegments), 3}; // [nSegments, 3] for x, y, z
                // H5::DataSpace posSpace(2, posDims);
                // H5::DataSet posDataset = liftingBladeGroup.createDataSet("pos_bound", H5::PredType::NATIVE_DOUBLE, posSpace);

                // // 分配存储位置数据的向量
                // std::vector<double> posData(turbineParams.nSegments * 3); // [nSegments * 3] for x, y, z
                // for (int i = 0; i < turbineParams.nSegments; ++i)
                // {
                //     Vec3 position_bound = pos.boundAt(b, timestep, i); // 获取位置 (x, y, z)
                //     posData[i * 3 + 0] = position_bound.x;                    // x 分量
                //     posData[i * 3 + 1] = position_bound.y;                    // y 分量
                //     posData[i * 3 + 2] = position_bound.z;                    // z 分量
                //     // 调试输出
                //     // std::cout << "pos_bound(b=" << b << ", t=" << timestep << ", i=" << i
                //     //           << "): x=" << position_bound.x << ", y=" << position_bound.y << ", z=" << position_bound.z << std::endl;
                // }
                // posDataset.write(&posData[0], H5::PredType::NATIVE_DOUBLE);
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
            // 打开或创建 HDF5 文件
            H5::H5File file(outputFile, H5F_ACC_RDWR);

            // 1. 创建 /config 组
            H5::Group configGroup = file.createGroup("/config");

            // 1. 写入几何数据
            H5::Group geomGroup = configGroup.createGroup("geometry");

            // 写入 rShedding [nSegments]
            hsize_t rDims[1] = {geom.rShedding.size()};
            H5::DataSpace rSpace(1, rDims);
            H5::DataSet rDataset = geomGroup.createDataSet("r_shed", H5::PredType::NATIVE_DOUBLE, rSpace);
            rDataset.write(geom.rShedding.data(), H5::PredType::NATIVE_DOUBLE);

            // 写入 rTrail [nSegments]
            hsize_t rTDims[1] = {geom.rTrailing.size()};
            H5::DataSpace rTSpace(1, rTDims);
            H5::DataSet rTDataset = geomGroup.createDataSet("r_trail", H5::PredType::NATIVE_DOUBLE, rTSpace);
            rTDataset.write(geom.rTrailing.data(), H5::PredType::NATIVE_DOUBLE);

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
            H5::Group simGroup = configGroup.createGroup("simulation");

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

    // 读取HDF5快照
    void read_wake_snapshot(Wake &wake, const std::string &h5_filepath, int timestep, const TurbineParams &turbineParams)
    {
        try
        {
            // 以只读模式打开HDF5文件
            H5::H5File file(h5_filepath, H5F_ACC_RDONLY);

            // --- 核心修复：在填充数据之前，确保Wake对象内部的vector有足够的空间 ---
            wake.ensureTimeStepExists(timestep);

            std::string wakeTimestepName = "/wake/timestep_" + std::to_string(timestep);
            if (!file.exists(wakeTimestepName))
            {
                throw std::runtime_error("Timestep " + std::to_string(timestep) + " not found in HDF5 file.");
            }

            for (int b = 0; b < turbineParams.nBlades; ++b)
            {
                BladeWake &bladeWake = wake.getBladeWake(timestep, b);
                std::string bladeGroupName = wakeTimestepName + "/blade_" + std::to_string(b);

                // 1. 读取节点数据
                std::string nodeDatasetName = bladeGroupName + "/nodes";
                if (file.exists(nodeDatasetName))
                {
                    H5::DataSet nodeDataset = file.openDataSet(nodeDatasetName);
                    H5::DataSpace nodeSpace = nodeDataset.getSpace();

                    hsize_t nodeDims[2];
                    nodeSpace.getSimpleExtentDims(nodeDims, NULL);
                    size_t nNodes = nodeDims[0];

                    std::vector<double> nodeData(nNodes * 6);
                    nodeDataset.read(nodeData.data(), H5::PredType::NATIVE_DOUBLE);

                    bladeWake.nodes.resize(nNodes);
                    for (size_t i = 0; i < nNodes; ++i)
                    {
                        bladeWake.nodes[i].position = {nodeData[i * 6 + 0], nodeData[i * 6 + 1], nodeData[i * 6 + 2]};
                        bladeWake.nodes[i].velocity = {nodeData[i * 6 + 3], nodeData[i * 6 + 4], nodeData[i * 6 + 5]};
                    }
                }

                // 2. 读取涡线数据
                std::string lineDatasetName = bladeGroupName + "/lines";
                if (file.exists(lineDatasetName))
                {
                    H5::DataSet lineDataset = file.openDataSet(lineDatasetName);
                    H5::DataSpace lineSpace = lineDataset.getSpace();

                    hsize_t lineDims[2];
                    lineSpace.getSimpleExtentDims(lineDims, NULL);
                    size_t nLines = lineDims[0];

                    std::vector<double> lineData(nLines * 4);
                    lineDataset.read(lineData.data(), H5::PredType::NATIVE_DOUBLE);

                    bladeWake.lines.resize(nLines);
                    for (size_t i = 0; i < nLines; ++i)
                    {
                        bladeWake.lines[i].startNodeIdx = static_cast<int>(lineData[i * 4 + 0]);
                        bladeWake.lines[i].endNodeIdx = static_cast<int>(lineData[i * 4 + 1]);
                        bladeWake.lines[i].gamma = lineData[i * 4 + 2];
                        bladeWake.lines[i].type = static_cast<VortexLineType>(lineData[i * 4 + 3]);
                    }
                }
            }
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "HDF5 read error: " << e.getDetailMsg() << std::endl;
        }
    }
    // 读取HDF5快照


    /**
     * @brief 将最终的拉格朗日尾流数据投影到欧拉网格上，并保存为VTK文件。
     */
    void projectWakeToEulerianGrid(const Wake &wake,
                                   const TurbineParams &turbineParams,
                                   const SimParams &simParams,
                                   const BladeGeometry &geom,
                                   int final_timestep,
                                   const std::string &outputPath)
    {
        // ===================== 参数设置 =====================
        std::filesystem::path path(outputPath);
        path /= "vel_eulerian.vtk";
        const std::string vtk_filepath = path.string();

        const double TURBINE_DIAMETER = turbineParams.rTip * 2.0;
        const double hub_height = turbineParams.hubHeight > 1.0 ? turbineParams.hubHeight : 90.0; // Use config or default

        // 网格参数选择
        const bool use_uniform_grid = true; // true: 均匀网格, false: 非均匀网格
        const double res_high_m = 1.0;       // 高分辨率区域的网格大小（米）

        // 网格范围设置
        const double x_start_m = -1.0 * TURBINE_DIAMETER;
        const double x_end_m = 7.0 * TURBINE_DIAMETER;
        const double y_max_m = 1.5 * TURBINE_DIAMETER / 2.0;
        const double z_max_m = 1.5 * TURBINE_DIAMETER / 2.0;

        // 非均匀网格参数 (如果 use_uniform_grid = false)
        const double res_low_m = 5.0;
        const double x_fine_start_m = 0.0 * TURBINE_DIAMETER;
        const double x_fine_end_m = 3.0 * TURBINE_DIAMETER;
        const double yz_fine_range_m = 1.5 * TURBINE_DIAMETER / 2.0;

        // ===================== 网格构建 =====================
        std::cout << "Creating Eulerian grid (" << (use_uniform_grid ? "uniform" : "non-uniform") << ")..." << std::endl;
        std::vector<double> x_coords, y_coords, z_coords;

        if (use_uniform_grid)
        {
            for (double x = x_start_m; x <= x_end_m; x += res_high_m)
                x_coords.push_back(x);
            for (double y = -y_max_m; y <= y_max_m; y += res_high_m)
                y_coords.push_back(y);
            for (double z = hub_height - z_max_m; z <= hub_height + z_max_m; z += res_high_m)
                z_coords.push_back(z);
        }
        else
        {
            // X coordinates
            for (double x = x_start_m; x < x_fine_start_m; x += res_low_m)
                x_coords.push_back(x);
            for (double x = x_fine_start_m; x < x_fine_end_m; x += res_high_m)
                x_coords.push_back(x);
            for (double x = x_fine_end_m; x <= x_end_m; x += res_low_m)
                x_coords.push_back(x);
            // Y coordinates
            for (double y = -y_max_m; y < -yz_fine_range_m; y += res_low_m)
                y_coords.push_back(y);
            for (double y = -yz_fine_range_m; y <= yz_fine_range_m; y += res_high_m)
                y_coords.push_back(y);
            for (double y = yz_fine_range_m + res_low_m; y <= y_max_m; y += res_low_m)
                y_coords.push_back(y);
            // Z coordinates
            for (double z = hub_height - z_max_m; z < hub_height - yz_fine_range_m; z += res_low_m)
                z_coords.push_back(z);
            for (double z = hub_height - yz_fine_range_m; z <= hub_height + yz_fine_range_m; z += res_high_m)
                z_coords.push_back(z);
            for (double z = hub_height + yz_fine_range_m + res_low_m; z <= hub_height + z_max_m; z += res_low_m)
                z_coords.push_back(z);
        }

        const int Nx = x_coords.size();
        const int Ny = y_coords.size();
        const int Nz = z_coords.size();

        std::vector<Vec3> grid_points;
        grid_points.reserve(Nx * Ny * Nz);
        for (int k = 0; k < Nz; ++k)
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    grid_points.push_back({x_coords[i], y_coords[j], z_coords[k]});

        std::cout << "Grid created. Dimensions: " << Nx << " x " << Ny << " x " << Nz
                  << ", Total points: " << grid_points.size() << std::endl;

        // ===================== 计算速度 =====================
        std::cout << "Calculating induced velocity on grid points..." << std::endl;
        std::vector<Vec3> grid_velocities;
        fvw::computeInducedVelocity(grid_velocities, grid_points, wake, final_timestep, turbineParams, geom, simParams);

        // ===================== 输出VTK =====================
        std::cout << "Writing velocity field to VTK file: " << vtk_filepath << std::endl;
        
        // Inline implementation for now to avoid header leakage
         std::ofstream vtkFile(vtk_filepath);
        if (!vtkFile.is_open())
        {
            std::cerr << "Error: Unable to open VTK file " << vtk_filepath << std::endl;
            return;
        }

        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "Wake Velocity Field\n";
        vtkFile << "ASCII\n";
        vtkFile << "DATASET STRUCTURED_GRID\n";
        vtkFile << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
        vtkFile << "POINTS " << grid_points.size() << " float\n";
        for (const auto &p : grid_points)
            vtkFile << p.x << " " << p.y << " " << p.z << "\n";

        vtkFile << "POINT_DATA " << grid_velocities.size() << "\n"; 
        vtkFile << "VECTORS Velocity float\n";
        for (const auto &v : grid_velocities) 
            vtkFile << v.x << " " << v.y << " " << v.z << "\n";

        vtkFile.close();
        
        std::cout << "Post-processing successfully completed!" << std::endl;
    }

    void runProbeCalculation(const std::string &h5_filepath, const std::string &csv_filepath,
                             const TurbineParams &turbineParams, const SimParams &simParams,
                             const BladeGeometry &geom)
    {
        double D = turbineParams.rTip * 2.0;

        // 2. 定义你的“虚拟探针”位置
        std::vector<Vec3> probe_points;
        double y_probe = 0.0;
        double z_probe = (turbineParams.hubHeight > 1.0 ? turbineParams.hubHeight : 90.0) + 50.0; // z=Hub+50m
        std::vector<double> x_locations_D = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
        for (double x_D : x_locations_D)
        {
            probe_points.push_back({x_D * D, y_probe, z_probe});
        }

        // 3. 准备输出文件
        std::ofstream outfile(csv_filepath);
        if (!outfile.is_open())
        {
            std::cerr << "Error: Unable to open CSV file for probe output: " << csv_filepath << std::endl;
            return;
        }
        outfile << "Timestep,Time,ProbeID,ProbeX,ProbeY,ProbeZ,U_ind,V_ind,W_ind\n";
        outfile << std::fixed << std::setprecision(6);

        // 4. 循环遍历HDF5中的所有时间步
        int start_step = 0;
        int end_step = simParams.timesteps - 1;
        int step_interval = simParams.outputFrequency;
        double dt = simParams.dt;

        // 创建一个Wake对象
        Wake wake(turbineParams.nBlades, turbineParams.nSegments, turbineParams.nSegments + 1);

        for (int t = start_step; t <= end_step; t += step_interval)
        {
            read_wake_snapshot(wake, h5_filepath, t, turbineParams);

            std::vector<Vec3> induced_velocities_at_probes;
            fvw::computeInducedVelocity(induced_velocities_at_probes, probe_points, wake, t, turbineParams, geom, simParams);

            for (size_t i = 0; i < probe_points.size(); ++i)
            {
                const auto &p = probe_points[i];
                const auto &vel = induced_velocities_at_probes[i];

                outfile << t << "," << t * dt << "," << i << ","
                        << p.x << "," << p.y << "," << p.z << ","
                        << vel.x << "," << vel.y << "," << vel.z << "\n";
            }
        }

        outfile.close();
        std::cout << "\n探针后处理完成，结果已保存至 " << csv_filepath << std::endl;
    }

} // namespace fvw