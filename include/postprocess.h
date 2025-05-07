#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "wake.h"
#include <string>

namespace fvw
{
    // 将 wake 数据写入 VTK 文件
    void writeWakeToVTK(const Wake &wake, const TurbineParams &turbineParams,
                        const std::string &outputDir, int timestep);
}

#endif // POSTPROCESS_H