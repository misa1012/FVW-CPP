#ifndef FVW_PERFORMANCE_H
#define FVW_PERFORMANCE_H

#include "utils.h"
#include <vector>
#include <stdexcept>
#include <iostream>

namespace fvw
{

    class PerformanceData
    {
    private:
        std::vector<double> aoa;              // 迎角 (deg)
        std::vector<double> cl;               // 升力系数
        std::vector<double> cd;               // 阻力系数
        std::vector<Vec3> inducedVelocity;    // 诱导速度 (叶片坐标系, u, v, w)
        std::vector<Vec3> inducedVelocityICS; // 诱导速度 (惯性坐标系)
        std::vector<double> boundGamma;       // Bound vortex strength (Gamma)
        int nBlades, nTimesteps, nShed;

    public:
        PerformanceData(int nBlades_, int timesteps_, int nShed_);

        // 只读访问器
        const double &aoaAt(int b, int t, int i) const;
        const double &clAt(int b, int t, int i) const;
        const double &cdAt(int b, int t, int i) const;
        const Vec3 &inducedVelocityAt(int b, int t, int i) const;
        const Vec3 &inducedVelocityICSAt(int b, int t, int i) const;

        // 可写访问器
        double &setAoaAt(int b, int t, int i);
        double &setClAt(int b, int t, int i);
        double &setCdAt(int b, int t, int i);
        Vec3 &setInducedVelocityAt(int b, int t, int i);
        Vec3 &setInducedVelocityICSAt(int b, int t, int i);

        // Getter 方法
        int getTimesteps() const { return nTimesteps; }
        int getBlades() const { return nBlades; }
        int getShed() const { return nShed; }

        // bound gamma
        double boundGammaAt(int b, int t, int i) const
        {
            return boundGamma[b * nTimesteps * nShed + t * nShed + i];
        }
        double &setBoundGammaAt(int b, int t, int i)
        {
            return boundGamma[b * nTimesteps * nShed + t * nShed + i];
        }
    };

} // namespace fvw

#endif // FVW_PERFORMANCE_H
