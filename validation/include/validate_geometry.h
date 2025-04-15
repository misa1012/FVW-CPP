#ifndef FVW_VALIDATE_GEOMETRY_H
#define FVW_VALIDATE_GEOMETRY_H

#include "geometry.h"

namespace fvw
{

    // 验证叶片几何并输出结果
    // saveToFile: 是否保存到 CSV
    void validateBladeGeometry(const BladeGeometry &geom, const TurbineParams &params, bool saveToFile = false);

} // namespace fvw

#endif // FVW_VALIDATE_GEOMETRY_H