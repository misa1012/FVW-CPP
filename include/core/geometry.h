#ifndef FVW_GEOMETRY_H
#define FVW_GEOMETRY_H

#include "core/utils.h"
#include "core/types.h"
#include <vector>

namespace fvw
{

    // Raw blade definition data (loaded from file)
    struct BladeDefinition
    {
        std::vector<double> r;
        std::vector<double> chord;
        std::vector<double> twist;
        std::vector<int> airfoilIndex;
    };

    struct BladeGeometry
    {
        std::vector<double> rTrailing;
        std::vector<double> chordTrailing;
        std::vector<double> twistTrailing;
        std::vector<double> rShedding;
        std::vector<double> chordShedding;
        std::vector<double> twistShedding;
        std::vector<Vec3> lead, quarter, trail; // Trailing nodes
        std::vector<Vec3> colloc, bound, end;   // Shedding nodes
        std::vector<int> airfoilIndex;
    };

    // 计算叶片几何
    // Now requires explicit BladeDefinition input instead of hardcoded internal data
    BladeGeometry computeBladeGeometry(const TurbineParams &params, const BladeDefinition &rawDist);

    // Reads blade geometry distribution from a CSV file
    // Format: r, chord, twist, airfoil_index
    BladeDefinition loadBladeDefinition(const std::string &csvPath);

} // namespace fvw

#endif // FVW_GEOMETRY_H
