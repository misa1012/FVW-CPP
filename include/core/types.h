#ifndef FVW_TYPES_H
#define FVW_TYPES_H

#include <string>
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
        VanGarrel,                // 0: 原始 Van Garrel (epsilon^2)
        VanGarrelUnitConsistent,  // 1: Van Garrel unit-consistent (epsilon^4)
        ChordBasedCore,           // 2: 跟local chord挂钩 建议选10%
    };

    // 定义节点分布类型
    enum class SegmentDistribution
    {
        Linear, // 线性均匀分布
        Cosine  // 余弦分布 (两端加密)
    };

    // 定义扰动类型
    enum class PerturbationType
    {
        None,                 // 无扰动
        CollectivePitch,      // 集体变桨 (所有叶片同步)
        AsymmetricStaticPitch // 非对称固定变桨
    };

    // 定义时间推进格式
    enum class TimeMarchingScheme
    {
        Euler,              // 一阶显式欧拉
        PredictorCorrector, // 二阶显式预测-校正 (Heun)
    };

    // Tip-loss correction model applied once at the sectional aerodynamic level.
    enum class TipLossModel
    {
        Off,
        Shen,
    };

    // *** 定义扰动参数 ***
    struct PerturbationParams
    {
        PerturbationType type;      // 默认无扰动
        double amplitude_deg = 0.0; // 扰动幅值 (度)
        double frequency_hz = 0.0;  // 扰动频率 (Hz)
    };

    // 仿真控制参数
    struct SimParams
    {
        double dt;
        double totalTime;
        int timesteps;
        int outputFrequency;
        int probeFrequency = 0;
        bool computeProbes = false;
        double cutoffParam;
        VortexModelType vortexModel;
        VortexCoreType coreType;
        PerturbationParams perturbation;
        
        // BEM Solver Settings
        double bemTolerance = 1e-4;
        int bemMaxIterations = 200;
        double bemRelaxation = 0.2;

        // Kutta-Joukowski iteration settings
        double kuttaTolerance = 1e-4;
        int kuttaMaxIterations = 200;
        double kuttaRelaxation = 0.3;
        
        // Revolution-based parameters (optional, for user convenience)
        // If > 0, dt and totalTime will be calculated automatically
        int stepsPerRevolution = 0;
        double numRevolutions = 0.0;

        // Logging controls
        bool logStepTiming = true;
        bool logVerbose = false;
        bool logPerf = false;

        // Time marching scheme for wake convection
        TimeMarchingScheme timeScheme = TimeMarchingScheme::Euler;

        // Single tip-loss correction switch (default off).
        TipLossModel tipLossModel = TipLossModel::Off;
        // Legacy coefficient fields retained for config compatibility; unused by
        // the current single-factor Shen implementation.
        double c1Faxi = 0.1219;
        double c2Faxi = 21.52;
        double c3Faxi = 0.1;
        double c1Ftan = 0.0984;
        double c2Ftan = 13.026;
        double c3Ftan = 0.1;
        
    };

    // 风机物理参数
    struct TurbineParams
    {
        std::string model;
        double windSpeed = -1.0;
        double rho = 1.225;
        double rHub = -1.0;
        double rTip = -1.0;
        double hubHeight = -1.0;
        int nBlades = -1;
        int nSegments = -1;
        SegmentDistribution segmentDistribution = SegmentDistribution::Linear;
        double tsr = -1.0;
        double omega = 0.0;
    };

} // namespace fvw

#endif // FVW_TYPES_H
