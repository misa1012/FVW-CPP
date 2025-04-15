#include "geometry.h"
#include "utils.h"
#include <cmath>

namespace fvw
{

    BladeGeometry computeBladeGeometry(const TurbineParams &params)
    {
        BladeGeometry geom;
        int nTrail = params.nSegments + 1;
        int nShed = params.nSegments;

        // Original blade data
        std::vector<double> originalR = {1.5, 2.86670, 5.6, 8.33330, 11.75, 15.85, 19.95,
                                         24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65,
                                         52.75, 56.1667, 58.90, 61.6333, 63.0};
        std::vector<double> originalChord = {3.542, 3.542, 3.854, 4.167, 4.557, 4.652, 4.458,
                                             4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764,
                                             2.518, 2.313, 2.086, 1.419, 1.419 / 2};
        std::vector<double> originalTwist = {-13.308, -13.308, -13.308, -13.308, -13.308, -11.48,
                                             -10.162, -9.011, -7.795, -6.544, -5.361, -4.188, -3.125,
                                             -2.319, -1.526, -0.863, -0.370, -0.106, -0.106 / 2};
        std::vector<int> originalNFoil = {0, 0, 0, 1, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7, 7};

        // Interpolation functions
        auto interpolateDouble = [](const std::vector<double> &x, const std::vector<double> &y, double xq)
        {
            if (xq <= x.front())
                return y.front();
            if (xq >= x.back())
                return y.back();
            for (size_t i = 1; i < x.size(); ++i)
            {
                if (xq <= x[i])
                {
                    double t = (xq - x[i - 1]) / (x[i] - x[i - 1]);
                    return y[i - 1] + t * (y[i] - y[i - 1]);
                }
            }
            return y.back();
        };

        auto interpolateInt = [](const std::vector<double> &x, const std::vector<int> &y, double xq)
        {
            if (xq <= x.front())
                return y.front();
            if (xq >= x.back())
                return y.back();
            for (size_t i = 1; i < x.size(); ++i)
            {
                if (xq <= x[i])
                {
                    // 比较 xq 到 x[i-1] 和 x[i] 的距离
                    double dist_prev = std::abs(xq - x[i - 1]);
                    double dist_curr = std::abs(xq - x[i]);
                    return (dist_prev <= dist_curr) ? y[i - 1] : y[i];
                }
            }
            return y.back();
        };

        // Cosine distribution for trailing nodes
        geom.rTrailing.resize(nTrail);
        geom.chordTrailing.resize(nTrail);
        geom.twistTrailing.resize(nTrail);
        std::vector<double> theta(nTrail);
        for (int i = 0; i < nTrail; ++i)
        {
            theta[i] = M_PI * i / (nTrail - 1);
            geom.rTrailing[i] = params.rHub + (params.rTip - params.rHub) * (1 - std::cos(theta[i])) / 2;
            geom.chordTrailing[i] = interpolateDouble(originalR, originalChord, geom.rTrailing[i]);
            geom.twistTrailing[i] = interpolateDouble(originalR, originalTwist, geom.rTrailing[i]);
        }

        // Shedding nodes
        geom.rShedding.resize(nShed);
        geom.chordShedding.resize(nShed);
        geom.twistShedding.resize(nShed);
        geom.airfoilIndex.resize(nShed);
        for (int i = 0; i < nShed; ++i)
        {
            geom.rShedding[i] = (geom.rTrailing[i] + geom.rTrailing[i + 1]) / 2;
            geom.chordShedding[i] = interpolateDouble(geom.rTrailing, geom.chordTrailing, geom.rShedding[i]);
            geom.twistShedding[i] = interpolateDouble(geom.rTrailing, geom.twistTrailing, geom.rShedding[i]);
            geom.airfoilIndex[i] = interpolateInt(originalR, originalNFoil, geom.rShedding[i]);
        }

        // Node positions
        // This is defined in BCS
        geom.lead.resize(nTrail);
        geom.quarter.resize(nTrail);
        geom.trail.resize(nTrail);
        for (int i = 0; i < nTrail; ++i)
        {
            geom.lead[i] = Vec3(-0.25 * geom.chordTrailing[i], 0.0, geom.rTrailing[i]);
            geom.quarter[i] = Vec3(0.0, 0.0, geom.rTrailing[i]);
            geom.trail[i] = Vec3(0.75 * geom.chordTrailing[i], 0.0, geom.rTrailing[i]);
        }

        geom.colloc.resize(nShed);
        geom.bound.resize(nShed);
        geom.end.resize(nShed);
        for (int i = 0; i < nShed; ++i)
        {
            geom.colloc[i] = Vec3(0.25 * geom.chordShedding[i], 0.0, geom.rShedding[i]);
            geom.bound[i] = Vec3(0.0, 0.0, geom.rShedding[i]);
            geom.end[i] = Vec3(0.75 * geom.chordShedding[i], 0.0, geom.rShedding[i]);
        }

        return geom;
    }

} // namespace fvw