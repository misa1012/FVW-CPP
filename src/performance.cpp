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

    // std::vector<double> computeAoA_tan(const PerformanceData &perf, const VelBCS &velBCS)
    // {
    //     std::vector<double> aoa(perf.getBlades() * perf.getTimesteps() * perf.getShed());
    //     for (int b = 0; b < perf.getBlades(); ++b)
    //     {
    //         for (int t = 0; t < perf.getTimesteps(); ++t)
    //         {
    //             for (int i = 0; i < perf.getShed(); ++i)
    //             {
    //                 int idx = b * perf.getTimesteps() * perf.getShed() + t * perf.getShed() + i;
    //                 aoa[idx] = std::atan2(-velBCS.at(b, t, i).y, velBCS.at(b, t, i).x) * 180.0 / M_PI;
    //                 if (std::isnan(aoa[idx]) || std::abs(aoa[idx]) == 180.0)
    //                 {
    //                     aoa[idx] = 0.0;
    //                 }
    //             }
    //         }
    //     }
    //     return aoa;
    // }

} // namespace fvw