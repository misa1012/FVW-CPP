#include "bem.h"
#include <algorithm>
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

    // 线性插值翼型数据
    static double interpolateAirfoil(const std::vector<double> &x,
                                     const std::vector<double> &y,
                                     double xq)
    {
        if (x.empty() || y.empty())
            return 0.0;
        if (xq <= x.front())
            return y.front();
        if (xq >= x.back())
            return y.back();

        for (size_t i = 1; i < x.size(); ++i)
        {
            if (xq <= x[i])
            {
                double t = (xq - x[i - 1]) / (x[i] - x[i - 1]);
                return y[i - 1] + t * (y[i] - y[i - 1]);
            }
        }
        return y.back();
    }

    void computeBEM(PerformanceData &perf,
                    const std::vector<std::vector<std::vector<double>>> &aoag,
                    const BladeGeometry &geom,
                    const TurbineParams &turbineParams,
                    const std::vector<AirfoilData> &airfoils)
    {
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int t = 0; t < perf.getTimesteps(); ++t)
            {
                for (int i = 0; i < turbineParams.nSegments; ++i)
                {
                    // 获取几何迎角
                    double aoa_deg = aoag[b][t][i];

                    // 限制迎角范围（与翼型数据匹配）
                    aoa_deg = std::max(-20.0, std::min(20.0, aoa_deg));

                    // 获取 shedding 节点的翼型索引
                    double rShed = geom.rShedding[i];
                    int af_idx = 0;
                    for (size_t j = 1; j < geom.rTrailing.size(); ++j)
                    {
                        if (rShed <= geom.rTrailing[j])
                        {
                            af_idx = geom.airfoilIndex[j - 1];
                            break;
                        }
                    }
                    if (rShed > geom.rTrailing.back())
                    {
                        af_idx = geom.airfoilIndex.back();
                    }

                    if (af_idx < 0 || af_idx >= static_cast<int>(airfoils.size()))
                    {
                        af_idx = 0; // 回退到 Cylinder1
                    }

                    // 插值升力和阻力系数
                    double cl = interpolateAirfoil(airfoils[af_idx].aoa,
                                                   airfoils[af_idx].cl,
                                                   aoa_deg);
                    double cd = interpolateAirfoil(airfoils[af_idx].aoa,
                                                   airfoils[af_idx].cd,
                                                   aoa_deg);

                    // 存储结果
                    perf.setAoaAt(b, t, i) = aoa_deg;
                    perf.setClAt(b, t, i) = cl;
                    perf.setCdAt(b, t, i) = cd;
                }
            }
        }
    }

} // namespace fvw