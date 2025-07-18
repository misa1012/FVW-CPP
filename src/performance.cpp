#include "performance.h"
#include <cmath>

namespace fvw
{

    PerformanceData::PerformanceData(int nBlades_, int timesteps_, int nShed_)
        : nBlades(nBlades_), nTimesteps(timesteps_), nShed(nShed_)
    {
        aoa.resize(nBlades * nTimesteps * nShed);
        cl.resize(nBlades * nTimesteps * nShed);
        cd.resize(nBlades * nTimesteps * nShed);
        inducedVelocity.resize(nBlades * nTimesteps * nShed, Vec3{0.0, 0.0, 0.0});
        inducedVelocityICS.resize(nBlades * nTimesteps * nShed, Vec3{0.0, 0.0, 0.0});
        relativeVelocity.resize(nBlades * nTimesteps * nShed, Vec3{0.0, 0.0, 0.0}); 
        boundGamma.resize(nBlades_ * nTimesteps * nShed_, 0.0);
    }

    const double &PerformanceData::aoaAt(int b, int t, int i) const
    {
        return aoa[b * nTimesteps * nShed + t * nShed + i];
    }

    const double &PerformanceData::clAt(int b, int t, int i) const
    {
        return cl[b * nTimesteps * nShed + t * nShed + i];
    }

    const double &PerformanceData::cdAt(int b, int t, int i) const
    {
        return cd[b * nTimesteps * nShed + t * nShed + i];
    }
    double &PerformanceData::setAoaAt(int b, int t, int i)
    {
        return aoa[b * nTimesteps * nShed + t * nShed + i];
    }

    double &PerformanceData::setClAt(int b, int t, int i)
    {
        return cl[b * nTimesteps * nShed + t * nShed + i];
    }

    double &PerformanceData::setCdAt(int b, int t, int i)
    {
        return cd[b * nTimesteps * nShed + t * nShed + i];
    }

    const Vec3 &PerformanceData::inducedVelocityAt(int b, int t, int i) const
    {
        return inducedVelocity[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PerformanceData::setInducedVelocityAt(int b, int t, int i)
    {
        return inducedVelocity[b * nTimesteps * nShed + t * nShed + i];
    }

    const Vec3 &PerformanceData::inducedVelocityICSAt(int b, int t, int i) const
    {
        return inducedVelocityICS[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PerformanceData::setInducedVelocityICSAt(int b, int t, int i)
    {
        return inducedVelocityICS[b * nTimesteps * nShed + t * nShed + i];
    }

    const Vec3 &PerformanceData::relativeVelocityAt(int b, int t, int i) const
    {
        return relativeVelocity[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PerformanceData::setRelativeVelocityAt(int b, int t, int i)
    {
        return relativeVelocity[b * nTimesteps * nShed + t * nShed + i];
    }

} // namespace fvw