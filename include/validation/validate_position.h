#ifndef FVW_VALIDATE_POSITION_H
#define FVW_VALIDATE_POSITION_H

#include "core/position.h"
#include "core/geometry.h"
#include "core/utils.h"

namespace fvw
{

    void validatePosition(const PositionData &pos, const SimParams &simParams,
                          const TurbineParams &turbineParams, const BladeGeometry &geom,
                          int blade = -1, int timestep = 0, const std::string &nodeType = "all",
                          bool saveToFile = false);

} // namespace fvw

#endif // FVW_VALIDATE_POSITION_H