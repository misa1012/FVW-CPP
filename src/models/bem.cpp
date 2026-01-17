#include "models/bem.h"
#include "io/logger.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <set>

namespace fvw
{
    void computeBEM(PerformanceData &perf,
                    const BladeGeometry &geom, const TurbineParams &turbineParams,
                    const std::vector<AirfoilData> &airfoils, const SimParams &simParams)
    {
        // Settings from SimParams
        const double tolBEM = simParams.bemTolerance;
        const int maxIterBEM = simParams.bemMaxIterations;
        const double weightFactor = simParams.bemRelaxation;

        bool if_InitialGuess = true;

        // 读取物理参数
        double windSpeed = turbineParams.windSpeed;
        double omega = turbineParams.omega;
        const double pi = M_PI;

        // 初始化局部变量
        std::vector<double> a(perf.getBlades() * perf.getTimesteps() * perf.getShed(), 0.0);
        std::vector<double> ap(perf.getBlades() * perf.getTimesteps() * perf.getShed(), 0.0);
        std::vector<double> a0(perf.getBlades() * perf.getTimesteps() * perf.getShed());
        std::vector<double> ap0(perf.getBlades() * perf.getTimesteps() * perf.getShed());
        std::vector<double> solidity(perf.getShed());

        // 计算 solidity：表示某一径向位置的局部叶片面积相对于环形面积的比例
        for (int i = 0; i < perf.getShed(); ++i)
        {
            solidity[i] = (turbineParams.nBlades * geom.chordShedding[i]) / (2 * pi * geom.rShedding[i]);
        }



        if (if_InitialGuess)
        {
            for (int b = 0; b < perf.getBlades(); ++b)
            {
                for (int t = 0; t < perf.getTimesteps(); ++t)
                {
                    for (int i = 0; i < perf.getShed(); ++i)
                    {
                        int idx = b * perf.getTimesteps() * perf.getShed() + t * perf.getShed() + i;
                        double lambda_r = turbineParams.omega * geom.rShedding[i] / turbineParams.windSpeed;
                        double twist = -geom.twistShedding[i] * pi / 180.0;
                        double expr = 4.0 - 4.0 * pi * lambda_r * solidity[i] +
                                      pi * lambda_r * lambda_r * solidity[i] * (8.0 * twist + pi * solidity[i]);
                        if (expr < 0.0 || std::isnan(expr))
                        {
                            a[idx] = 0.33;
                            std::cerr << "[Warning] Invalid sqrt_expr at i=" << i << ", expr=" << expr << ", set a=0.33" << std::endl;
                        }
                        else
                        {
                            a[idx] = 0.25 * (2.0 + pi * lambda_r * solidity[i] - std::sqrt(expr));
                        }
                        ap[idx] = 0.0;
                        if (std::isnan(a[idx]) || std::isinf(a[idx]))
                        {
                            a[idx] = 0.33;
                            std::cerr << "[Warning] Initial a NaN/Inf at i=" << i << ", set a=0.33" << std::endl;
                        }
                    }
                }
            }
        }

        // BEM 迭代
        for (int iter = 0; iter < maxIterBEM; ++iter)
        {
            bool converged = true;

            for (int b = 0; b < perf.getBlades(); ++b)
            {
                for (int t = 0; t < perf.getTimesteps(); ++t)
                {
                    for (int i = 0; i < perf.getShed(); ++i)
                    {
                        int idx = b * perf.getTimesteps() * perf.getShed() + t * perf.getShed() + i;

                        // 保存旧值
                        a0[idx] = a[idx];
                        ap0[idx] = ap[idx];

                        double r = geom.rShedding[i];
                        double twist = -1 * geom.twistShedding[i] * M_PI / 180.0; // 转换为弧度

                        if (std::isnan(a[idx]) || std::isnan(ap[idx]) || std::isinf(a[idx]) || std::isinf(ap[idx]))
                        {
                            a[idx] = 0.33;
                            ap[idx] = 0.0;
                            std::cerr << "[Warning] Reset a, ap at i=" << i << " due to NaN/Inf" << std::endl;
                        }

                        double phi = std::atan2(windSpeed * (1.0 - a[idx]), omega * r * (1.0 + ap[idx]));
                        // 计算迎角
                        double aoa = (phi - twist) * 180.0 / M_PI;

                        if (std::isnan(aoa) || std::abs(aoa) == 180.0)
                        {
                            aoa = 0.0;
                            std::cerr << "[Warning] AOA reset to 0 at i=" << i << ", phi=" << (phi * 180.0 / M_PI)
                                      << ", a=" << a[idx] << ", ap=" << ap[idx] << std::endl;
                        }
                        perf.setAoaAt(b, t, i) = aoa;



                        int airfoilIdx = geom.airfoilIndex[i];
                        perf.setClAt(b, t, i) = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cl, aoa);
                        perf.setCdAt(b, t, i) = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cd, aoa);

                        double cl = perf.clAt(b, t, i);
                        double cd = perf.cdAt(b, t, i);

                        // tip loss factor
                        double ftip = 2.0 / pi * std::acos(std::exp(-turbineParams.nBlades * (turbineParams.rTip - geom.rShedding[i]) / (2.0 * geom.rShedding[i] * std::sin(phi))));
                        double fhub = 2.0 / pi * std::acos(std::exp(-turbineParams.nBlades * (geom.rShedding[i] - turbineParams.rHub) / (2.0 * turbineParams.rHub * std::sin(phi))));
                        double f = ftip * fhub;
                        if (std::isnan(f) || f <= 0.0)
                        {
                            f = 1.0;
                            std::cerr << "[Warning] Reset f at i=" << i << " due to NaN/negative" << std::endl;
                        }

                        // 推力系数
                        double ct = solidity[i] * (1.0 - a[idx]) * (1.0 - a[idx]) *
                                    (cl * std::cos(phi) + cd * std::sin(phi)) /
                                    (std::sin(phi) * std::sin(phi));

                        // 更新 a
                        double denom = solidity[i] * (cl * std::cos(phi) + cd * std::sin(phi));
                        double a_new = denom > 1e-10 ? 1.0 / (1.0 + 4.0 * f * std::sin(phi) * std::sin(phi) / denom) : a[idx];
                        if (ct > 0.96 * f)
                        {
                            a_new = (18.0 * f - 20.0 - 3.0 * std::sqrt(ct * (50.0 - 36.0 * f) + 12.0 * f * (3.0 * f - 4.0))) /
                                    (36.0 * f - 50.0);
                        }

                        // 更新 ap
                        double ap_new = denom > 1e-10 ? 1.0 / (4.0 * f * std::cos(phi) * std::sin(phi) /
                                                                   (solidity[i] * (cl * std::sin(phi) - cd * std::cos(phi))) -
                                                               1.0)
                                                      : ap[idx];

                        // 残差
                        double da = std::abs(a_new - a0[idx]);
                        double dap = std::abs(ap_new - ap0[idx]);

                        // 加权更新
                        a[idx] = (1.0 - weightFactor) * a[idx] + weightFactor * a_new;
                        ap[idx] = (1.0 - weightFactor) * ap[idx] + weightFactor * ap_new;

                        if (da > tolBEM || dap > tolBEM)
                        {
                            converged = false;
                        }
                    }
                }
            }

            if (converged)
            {
                Logger::log("SOLVER", "BEM Converged: " + std::to_string(iter + 1) + " iters");

                // 所有迭代完毕了输出


                break;
            }
            if (iter == maxIterBEM - 1)
            {
                std::cerr << "Warning: BEM failed to converge after " << maxIterBEM << " iterations" << std::endl;
            }
        }
    }

} // namespace fvw
