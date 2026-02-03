#include "core/utils.h"

namespace fvw
{
    // 线性插值函数
    double interpolate(const std::vector<double> &x, const std::vector<double> &y, double x0)
    {
        if (x.empty() || x.size() != y.size())
            return 0.0;
        if (x0 <= x.front())
            return y.front();
        if (x0 >= x.back())
            return y.back();

        for (size_t i = 0; i < x.size() - 1; ++i)
        {
            if (x[i] <= x0 && x0 <= x[i + 1])
            {
                double t = (x0 - x[i]) / (x[i + 1] - x[i]);
                return y[i] + t * (y[i + 1] - y[i]);
            }
        }
        return 0.0;
    }

    // 用于airfoil插值
    int interpolateInt(const std::vector<double> &x, const std::vector<int> &y, double xq)
    {
        // 1. 输入验证
        if (x.empty() || x.size() != y.size())
        {
            return 0; // 对于 int 类型，返回 0 作为错误/无效输入的指示
        }

        // 2. 边界检查
        if (xq <= x.front())
        {
            return y.front();
        }
        if (xq >= x.back())
        {
            return y.back();
        }

        // 3. 循环查找区间并进行最近邻插值
        for (size_t i = 0; i < x.size() - 1; ++i)
        {
            // 检查 xq 是否在当前区间 [x[i], x[i+1]] 内
            // (假设 x 是已排序的)
            if (x[i] <= xq && xq <= x[i + 1])
            {
                double dist_to_xi = std::abs(xq - x[i]);
                double dist_to_xi_plus_1 = std::abs(xq - x[i + 1]);

                // 如果 xq 更接近或等距于 x[i]，则返回 y[i]
                // 否则返回 y[i+1]
                if (dist_to_xi <= dist_to_xi_plus_1)
                {
                    return y[i];
                }
                else
                {
                    return y[i + 1];
                }
            }
        }

        // 4. 回退返回值
        //    理论上，如果 x 有序且 xq 在 x 的范围内，不应该执行到这里。
        //    返回 0 作为 int 类型的回退值。
        return 0;
    }

    std::string to_string(const Vec3 &v)
    {
        return "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")";
    }

} // namespace fvw