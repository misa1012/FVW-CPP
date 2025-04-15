#ifndef FVW_VALIDATE_H
#define FVW_VALIDATE_H

#include "geometry.h"
#include "position.h"
#include "utils.h"
#include <string>
#include <vector>

namespace fvw
{

    struct ValidationParams
    {
        std::string section = "all"; // "geometry", "position", "all"
        int blade = -1;
        int timestep = 0;
        std::string nodeType = "all";
        bool saveToFile = false;
    };

    ValidationParams parseValidationArgs(int argc, char *argv[]);

    void validate(const ValidationParams &params, const SimParams &simParams,
                  const TurbineParams &turbineParams, const BladeGeometry &geom,
                  const PositionData &pos);

} // namespace fvw

#endif // FVW_VALIDATE_H