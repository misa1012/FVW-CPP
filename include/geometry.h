#ifndef FVW_GEOMETRY_H
#define FVW_GEOMETRY_H

#include "utils.h"
#include <vector>

namespace fvw
{

    struct SimParams
    {
        double dt;
        double totalTime;
        int timesteps;
    };

    struct TurbineParams
    {
        double windSpeed;
        double rho;
        double rHub;
        double rTip;
        int nBlades;
        int nSegments;
        double tsr;
        double omega;
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
    };

    // 计算叶片几何
    BladeGeometry computeBladeGeometry(const TurbineParams &params);

} // namespace fvw

#endif // FVW_GEOMETRY_H