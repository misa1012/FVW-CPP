#ifndef FVW_PERFORMANCE_H
#define FVW_PERFORMANCE_H


#include "utils.h"
#include <vector>

namespace fvw
{

    class PerformanceData
    {
    private:
        std::vector<double> aoa; // 迎角 (deg)
        std::vector<double> cl;  // 升力系数
        std::vector<double> cd;  // 阻力系数
        int nBlades, nTimesteps, nShed;

    public:
        PerformanceData(int nBlades_, int timesteps_, int nShed_);

        // 只读访问器
        const double &aoaAt(int b, int t, int i) const;
        const double &clAt(int b, int t, int i) const;
        const double &cdAt(int b, int t, int i) const;

        // 可写访问器
        double &setAoaAt(int b, int t, int i);
        double &setClAt(int b, int t, int i);
        double &setCdAt(int b, int t, int i);

        // Getter 方法
        int getTimesteps() const { return nTimesteps; }
        int getBlades() const { return nBlades; }
        int getShed() const { return nShed; }
    };

    // std::vector<double> computeAoA_tan(const PerformanceData &perf, const VelBCS &velBCS);

} // namespace fvw

#endif // FVW_PERFORMANCE_H