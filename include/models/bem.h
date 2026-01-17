#ifndef FVW_BEM_H
#define FVW_BEM_H

#include "models/performance.h"
#include "core/geometry.h"
#include "models/airfoil.h"
#include "core/utils.h"

namespace fvw
{
    void computeBEM(PerformanceData &perf,
                    const BladeGeometry &geom, const TurbineParams &turbineParams,
                    const std::vector<AirfoilData> &airfoils);

} // namespace fvw

#endif // FVW_BEM_H