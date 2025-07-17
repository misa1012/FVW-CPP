#ifndef FVW_GEOMETRY_H
#define FVW_GEOMETRY_H

#include "utils.h"
#include <vector>

namespace fvw
{
    // 定义涡量模型的类型
    enum class VortexModelType
    {
        Constant,   // 0: 原始的、保持涡量不变
        GammaDecay, // 1: Gamma人工耗散+移除模型
    };

    // 定义涡核
    enum class VortexCoreType
    {
        VanGarrel,      // 0: 原始的、基于l0的平滑模型
        ChordBasedCore, // 1: 跟local chord挂钩 建议选10%
    };

    // 定义扰动类型
    enum class PerturbationType
    {
        None,                 // 无扰动
        CollectivePitch,      // 集体变桨 (所有叶片同步)
        AsymmetricStaticPitch // 非对称固定变桨
    };

    // *** 定义扰动参数 ***
    struct PerturbationParams
    {
        PerturbationType type;      // 默认无扰动
        double amplitude_deg = 0.0; // 扰动幅值 (度)
        double frequency_hz = 0.0;  // 扰动频率 (Hz)
    };

    struct SimParams
    {
        double dt;
        double totalTime;
        int timesteps;
        int outputFrequency;
        double cutoffParam;
        VortexModelType vortexModel;
        VortexCoreType coreType;
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
