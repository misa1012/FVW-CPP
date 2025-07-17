#ifndef FVW_WAKE_H
#define FVW_WAKE_H

#include "geometry.h"
#include "performance.h"
#include "position.h"
#include "utils.h"
#include "velocity.h"
#include "airfoil.h"
#include <vector>
#include <stdexcept>
#include <string> // 用于 to_string 等

namespace fvw
{

    struct VortexNode
    {
        Vec3 position;
        Vec3 velocity;
    };

    enum class VortexLineType
    {
        Bound,
        Shed,
        Trailing
    };

    // 需要包含vortex model相关的信息
    struct VortexLine
    {
        int startNodeIdx;
        int endNodeIdx;
        double gamma;
        VortexLineType type;

        // --- For Gamma Decay Model ---
        double initial_gamma; // 涡元诞生时的初始环量
        bool in_far_wake;     // 是否已进入远场开始衰减的标志

        // For chord base vortex core
        int segment_index;
    };

    // --- 单个叶片在一个时间步的尾迹信息 ---
    struct BladeWake
    {
        std::vector<VortexNode> nodes; // 该叶片在该时间步的所有尾迹节点 (包括叶片上的附着涡节点)
        std::vector<VortexLine> lines; // 连接该叶片内部节点的涡线段

        std::vector<int> boundNodeIndices; // t=n 时刻的附着涡节点索引 (nTrail)
        std::vector<int> trailNodeIndices; // t=n 时刻的尾迹节点索引 (nTrail)

        std::vector<double> prevGammaBound; // 存储 t=n-1 的 Bound 涡量强度

        // 存储不同类型涡线在 lines 向量中的索引，方便快速访问
        std::vector<int> boundLineIndices;    // Size nShed
        std::vector<int> trailingLineIndices; // Size nTrail
        std::vector<int> shedLineIndices;     // Size nShed

        // 构造函数
        BladeWake(int nTrail)
        {
            boundNodeIndices.resize(nTrail, -1);
            trailNodeIndices.resize(nTrail, -1);
            nodes.reserve(3 * nTrail);                    // 2*nTrail (对流) + nTrail (新附着涡)
            lines.reserve(3 * (nTrail - 1) + 2 * nTrail); // Bound + Trailing + Shed
            prevGammaBound.resize(nTrail - 1, 0.0);       // nShed = nTrail - 1

            // 初始化新增的索引向量
            boundLineIndices.resize(nTrail - 1, -1);
            trailingLineIndices.resize(nTrail, -1);
            shedLineIndices.resize(nTrail - 1, -1);
        }

        // 可以在这里添加辅助函数，例如添加节点/线段并返回索引
        int addNode(const VortexNode &node)
        {
            nodes.push_back(node);
            return nodes.size() - 1; // 返回新添加节点的索引
        }

        int addLine(const VortexLine &line)
        {
            // 可以添加检查，确保索引有效
            if (line.startNodeIdx >= nodes.size() || line.endNodeIdx >= nodes.size() || line.startNodeIdx < 0 || line.endNodeIdx < 0)
            {
                throw std::out_of_range("VortexLine node index out of range for BladeWake.");
            }
            lines.push_back(line);
            return lines.size() - 1; // 返回新添加线段的索引
        }
    };

    struct Wake
    {
        // 核心数据结构: [时间步][叶片] -> BladeWake
        std::vector<std::vector<BladeWake>> bladeWakes;

        // 参数信息 (方便访问)
        int nBlades;
        int nShed;
        int nTrail;

        // 构造函数
        Wake(int nBlades_, int nShed_, int nTrail_)
            : nBlades(nBlades_), nShed(nShed_), nTrail(nTrail_) {}

        // --- 辅助函数 ---
        // 确保指定时间步的结构存在，如果不存在则创建
        void ensureTimeStepExists(int timestep)
        {
            if (timestep >= bladeWakes.size())
            {
                bladeWakes.resize(timestep + 1); // 扩展外层 vector
            }
            if (bladeWakes[timestep].empty())
            {
                bladeWakes[timestep].resize(nBlades, BladeWake(nTrail)); // 为该时间步创建 nBlades 个 BladeWake
            }
        }

        // 获取特定叶片在特定时间步的尾迹数据 (非 const 版本)s
        BladeWake &getBladeWake(int timestep, int blade)
        {
            if (timestep >= bladeWakes.size() || blade < 0 || blade >= nBlades)
            {
                throw std::out_of_range("Timestep or blade index out of range in getBladeWake.");
            }
            if (bladeWakes[timestep].empty())
            {
                // 或者调用 ensureTimeStepExists(timestep);
                throw std::runtime_error("Accessing uninitialized timestep structure in getBladeWake.");
            }
            return bladeWakes.at(timestep).at(blade);
        }

        // 获取特定叶片在特定时间步的尾迹数据 (const 版本)
        const BladeWake &getBladeWake(int timestep, int blade) const
        {
            if (timestep >= bladeWakes.size() || blade < 0 || blade >= nBlades)
            {
                throw std::out_of_range("Timestep or blade index out of range in getBladeWake (const).");
            }
            if (bladeWakes[timestep].empty())
            {
                throw std::runtime_error("Accessing uninitialized timestep structure in getBladeWake (const).");
            }
            return bladeWakes.at(timestep).at(blade);
        }

        // 获取节点 (非 const)
        std::vector<VortexNode> &getNodes(int timestep, int blade)
        {
            return getBladeWake(timestep, blade).nodes;
        }

        // 获取节点 (const)
        const std::vector<VortexNode> &getNodes(int timestep, int blade) const
        {
            return getBladeWake(timestep, blade).nodes;
        }

        // 获取线段 (非 const)
        std::vector<VortexLine> &getLines(int timestep, int blade)
        {
            return getBladeWake(timestep, blade).lines;
        }
        // 获取线段 (const)
        const std::vector<VortexLine> &getLines(int timestep, int blade) const
        {
            return getBladeWake(timestep, blade).lines;
        }
    };

    // --- 函数声明 ---

    // Update the vortex core
    // void UpdateVortexStates(Wake &wake, int timestep, const SimParams &simParams, const TurbineParams &turbineParams);

    // Biot-Savart function
    void computeInducedVelocity(std::vector<Vec3> &inducedVelocities, const std::vector<Vec3> &targetPoints,
                                const Wake &wake, int timestep, const TurbineParams &turbineParams, const BladeGeometry &geom, const SimParams &simParams);

    void InitializeWakeStructure(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                                 const TurbineParams &turbineParams, const PositionData &pos, const SimParams &simParams);

    void UpdateWakeVelocities(Wake &wake, const TurbineParams &turbineParams, int timestep, const BladeGeometry &geom, const SimParams &simParams);

    void AdvanceWakeStructure(Wake &wake, const BladeGeometry &geom,
                              const TurbineParams &turbineParams, const PositionData &pos, double dt, int currentTimestep);

    void kuttaJoukowskiIteration(Wake &wake, PerformanceData &perf, const BladeGeometry &geom, NodeAxes &axes,
                                 const TurbineParams &turbineParams, const PositionData &pos, VelBCS &velBCS,
                                 std::vector<AirfoilData> &airfoils, const SimParams &simParams);

    void ApplyGammaDecayAndRemoval(Wake &wake, int timestep, const TurbineParams &turbineParams, const SimParams &simParams);
} // namespace fvw

#endif // FVW_WAKE_H