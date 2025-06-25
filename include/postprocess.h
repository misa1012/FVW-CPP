#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "wake.h"
#include "position.h"
#include <string>
#include <H5Cpp.h>

namespace fvw
{
    // 将 wake 数据写入 VTK 文件
    void writeWakeToVTK(const Wake &wake, const TurbineParams &turbineParams,
                        const std::string &outputDir, int timestep);

    // 把wake相关的信息写入hdf5
    void writeWakeToHDF5(const Wake &wake, const PositionData &pos, const PerformanceData &perf, VelICS &velICS, VelBCS &velBCS,
                         const TurbineParams &turbineParams, const std::string &outputFile,
                         int timestep);

    void writeConfigToHDF5(const BladeGeometry &geom, const TurbineParams &turbineParams,
                           const SimParams &simParams, const std::string &outputFile);

    void read_wake_snapshot(Wake &wake, const std::string &h5_filepath, int timestep, const TurbineParams &turbineParams);
}

#endif // POSTPROCESS_H