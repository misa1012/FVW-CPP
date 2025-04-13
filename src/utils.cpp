#include "utils.h"
#include <stdexcept>

namespace fvw
{

    // TODO: This is for one value; need to interplote for array
    double interpolate(const std::vector<double> &x, const std::vector<double> &y, double x_new)
    {
        if (x.size() != y.size() || x.empty())
        {
            throw std::invalid_argument("Invalid input for interpolation");
        }
        for (size_t i = 0; i < x.size() - 1; ++i)
        {
            if (x_new >= x[i] && x_new <= x[i + 1])
            {
                double t = (x_new - x[i]) / (x[i + 1] - x[i]);
                return y[i] + t * (y[i + 1] - y[i]);
            }
        }
        return x_new < x[0] ? y[0] : y.back();
    }

} // namespace fvw