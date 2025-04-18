#ifndef FVW_BEM_H
#define FVW_BEM_H

#include "performance.h"
#include "geometry.h"
#include "airfoil.h"
#include "utils.h"

namespace fvw
{

    void computeAoA(PerformanceData &perf,
                    const BladeGeometry &geom, const TurbineParams &turbineParams,
                    const std::vector<double> &a, const std::vector<double> &ap);

    void computeBEM(PerformanceData &perf,
                    const BladeGeometry &geom, const TurbineParams &turbineParams,
                    const std::vector<AirfoilData> &airfoils);

} // namespace fvw

#endif // FVW_BEM_H