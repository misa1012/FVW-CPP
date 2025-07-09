#ifndef FVW_GEOMETRY_H
#define FVW_GEOMETRY_H

#include "utils.h"
#include <vector>

namespace fvw
{
    // 定义涡核模型的类型
    enum class VortexModelType {
        VanGarrel,              // 0: 原始的、基于l0的平滑模型
        GammaDecay,             // 1: Gamma人工耗散+移除模型
        // QBlade                  // 2: QBlade的、带物理扩散和拉伸的高级模型
    };

    // 定义扰动类型
    enum class PerturbationType {
        None,               // 无扰动
        CollectivePitch     // 集体变桨 (所有叶片同步)
        // 未来可以扩展: CyclicPitch
    };

    // *** 定义扰动参数 ***
    struct PerturbationParams {
        PerturbationType type; // 默认无扰动
        double amplitude_deg = 0.0; // 扰动幅值 (度)
        double frequency_hz = 0.0;  // 扰动频率 (Hz)
    };

    struct SimParams
    {
        double dt;
        double totalTime;
        int timesteps;
        int outputFrequency;
        VortexModelType vortexModel;
        PerturbationParams perturbation;
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
        std::vector<int> airfoilIndex;
    };

    // 计算叶片几何
    BladeGeometry computeBladeGeometry(const TurbineParams &params);

} // namespace fvw

#endif // FVW_GEOMETRY_H
