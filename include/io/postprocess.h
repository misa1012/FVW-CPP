#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "core/wake.h"
#include "core/position.h"
#include <string>
#include <H5Cpp.h>

namespace fvw
{
    // 将 wake 数据写入 VTK 文件
    void writeWakeToVTK(const Wake &wake,
                        const std::string &outputDir, int timestep);

    // 把wake相关的信息写入hdf5
    void writeWakeToHDF5(const Wake &wake, const PerformanceData &perf,
                         const TurbineParams &turbineParams, const std::string &outputFile,
                         int timestep);

    void writeConfigToHDF5(const BladeGeometry &geom, const TurbineParams &turbineParams,
                           const SimParams &simParams, const std::string &outputFile);

    // Legacy/Main Integration functions (Moved from main.cpp)
    void projectWakeToEulerianGrid(const Wake &wake,
                                   const TurbineParams &turbineParams,
                                   const SimParams &simParams,
                                   const BladeGeometry &geom,
                                   int final_timestep,
                                   const std::string &outputPath);

    void runProbeCalculation(const std::string &h5_filepath, const std::string &csv_filepath,
                             const TurbineParams &turbineParams, const SimParams &simParams,
                             const BladeGeometry &geom);

    void read_wake_snapshot(Wake &wake, const std::string &h5_filepath, int timestep, const TurbineParams &turbineParams);
}

#endif // POSTPROCESS_H