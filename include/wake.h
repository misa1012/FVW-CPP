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
        int idx;
    };

    // 这个区分似乎，不太合理，区分了NewShed，也应该区分NewTrailing？
    enum class VortexLineType
    {
        Bound,
        Shed,
        Trailing
    };

    struct VortexLine
    {
        int startNodeIdx;
        int endNodeIdx;
        double gamma;
        VortexLineType type;
    };

    // --- 单个叶片在一个时间步的尾迹信息 ---
    struct BladeWake
    {
        std::vector<VortexNode> nodes;     // 该叶片在该时间步的所有尾迹节点 (包括叶片上的附着涡节点)
        std::vector<VortexLine> lines;     // 连接该叶片内部节点的涡线段
        std::vector<int> boundNodeIndices; // t=n 时刻的附着涡节点索引 (nTrail)
        std::vector<int> trailNodeIndices; // t=n 时刻的尾迹节点索引 (nTrail)
        std::vector<double> prevGammaBound; // 存储 t=n-1 的 Bound 涡量强度

        // 构造函数
        BladeWake(int nTrail)
        {
            boundNodeIndices.resize(nTrail, -1);
            trailNodeIndices.resize(nTrail, -1);
            nodes.reserve(3 * nTrail);                    // 2*nTrail (对流) + nTrail (新附着涡)
            lines.reserve(3 * (nTrail - 1) + 2 * nTrail); // Bound + Trailing + Shed
            prevGammaBound.resize(nTrail - 1, 0.0); // nShed = nTrail - 1
        }

        // 可以在这里添加辅助函数，例如添加节点/线段并返回索引
        int addNode(const VortexNode &node)
        {
            nodes.push_back(node);
            return nodes.size() - 1; // 返回新添加节点的索引
        }

        void addLine(const VortexLine &line)
        {
            // 可以添加检查，确保索引有效
            if (line.startNodeIdx >= nodes.size() || line.endNodeIdx >= nodes.size() || line.startNodeIdx < 0 || line.endNodeIdx < 0)
            {
                throw std::out_of_range("VortexLine node index out of range for BladeWake.");
            }
            lines.push_back(line);
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

        // 添加一个新的空时间步结构 (用于下一个时间步的初始化)
        // void addTimeStepStructure()
        // {
        //     bladeWakes.emplace_back(nBlades); // 添加 nBlades 个空的 BladeWake
        // }

        // 获取特定叶片在特定时间步的尾迹数据 (非 const 版本)
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

    void computeInducedVelocity(std::vector<Vec3> &inducedVelocities, const std::vector<Vec3> &targetPoints,
                                const Wake &wake, int timestep, const TurbineParams &turbineParams, double cutOff = 0.001);

    void InitializeWake(Wake &wake, const BladeGeometry &geom, PerformanceData &perf,
                        const TurbineParams &turbineParams, const PositionData &pos, double dt);

    void UpdateWakeVelocities(Wake &wake, const TurbineParams &turbineParams, int timestep);

    // void kuttaJoukowskiIteration(Wake &wake, PerformanceData &perf, const BladeGeometry &geom, NodeAxes &axes,
    //                              const TurbineParams &turbineParams, const PositionData &pos, VelBCS &velBCS, std::vector<AirfoilData> &airfoils);

} // namespace fvw

#endif // FVW_WAKE_H