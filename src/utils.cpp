#include "utils.h"

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

    std::string to_string(const Vec3& v) {
        return "(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")";
    }

} // namespace fvw