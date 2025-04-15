#ifndef FVW_BEM_H
#define FVW_BEM_H

#include "utils.h"
#include "geometry.h"
#include "airfoil.h"
#include <vector>

namespace fvw {

class PerformanceData {
private:
    std::vector<double> aoa; // 迎角 (deg)
    std::vector<double> cl;  // 升力系数
    std::vector<double> cd;  // 阻力系数
    int nBlades, nTimesteps, nShed;

public:
    PerformanceData(int nBlades_, int timesteps_, int nShed_);
    
    // 只读访问器
    const double& aoaAt(int b, int t, int i) const;
    const double& clAt(int b, int t, int i) const;
    const double& cdAt(int b, int t, int i) const;
    
    // 可写访问器
    double& setAoaAt(int b, int t, int i);
    double& setClAt(int b, int t, int i);
    double& setCdAt(int b, int t, int i);
    
    // Getter 方法
    int getTimesteps() const { return nTimesteps; }
    int getBlades() const { return nBlades; }
    int getShed() const { return nShed; }
};

// 计算 BEM 性能
void computeBEM(PerformanceData& perf,
                const std::vector<std::vector<std::vector<double>>>& aoag,
                const BladeGeometry& geom,
                const TurbineParams& turbineParams,
                const std::vector<AirfoilData>& airfoils);

} // namespace fvw

#endif // FVW_BEM_H
