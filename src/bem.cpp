#include "bem.h"
#include <iostream>

namespace fvw
{
    void computeBEM(PerformanceData &perf,
                    const BladeGeometry &geom, const TurbineParams &turbineParams,
                    const std::vector<AirfoilData> &airfoils)
    {
        // 定义一些开关
        bool if_InitialGuess = false;

        const double tolBEM = 1e-4;
        const int maxIterBEM = 200;
        const double pi = M_PI;
        const double weightFactor = 0.2;

        // 初始化局部变量
        std::vector<double> a(perf.getBlades() * perf.getTimesteps() * perf.getShed(), 0.0);
        std::vector<double> ap(perf.getBlades() * perf.getTimesteps() * perf.getShed(), 0.0);
        std::vector<double> a0(perf.getBlades() * perf.getTimesteps() * perf.getShed());
        std::vector<double> ap0(perf.getBlades() * perf.getTimesteps() * perf.getShed());
        std::vector<double> solidity(perf.getShed());
        std::vector<double> rNodes(perf.getShed());

        // 计算 solidity, rNodes
        for (int i = 0; i < perf.getShed(); ++i)
        {
            solidity[i] = (turbineParams.nBlades * geom.chordShedding[i]) / (2 * pi * geom.rShedding[i]);
            rNodes[i] = geom.rShedding[i];
        }

        // BEM 迭代
        for (int iter = 0; iter < maxIterBEM; ++iter)
        {
            bool converged = true;

            // 保存旧值
            for (int b = 0; b < perf.getBlades(); ++b)
            {
                for (int t = 0; t < perf.getTimesteps(); ++t)
                {
                    for (int i = 0; i < perf.getShed(); ++i)
                    {
                        int idx = b * perf.getTimesteps() * perf.getShed() + t * perf.getShed() + i;

                        if (if_InitialGuess)
                        {
                            double lambda_r = turbineParams.omega * rNodes[i] / turbineParams.windSpeed;
                            double twist = geom.twistShedding[i] * pi / 180.0;
                            a[idx] = 0.25 * (2.0 + pi * lambda_r * solidity[i] -
                                             std::sqrt(4.0 - 4.0 * pi * lambda_r * solidity[i] +
                                                       pi * lambda_r * lambda_r * solidity[i] *
                                                           (8.0 * twist + pi * solidity[i])));
                        }

                        a0[idx] = a[idx];
                        ap0[idx] = ap[idx];

                        double windSpeed = turbineParams.windSpeed;
                        double omega = turbineParams.omega;
                        double r = geom.rShedding[i];
                        double twist = geom.twistShedding[i] * M_PI / 180.0; // 转换为弧度

                        // 计算流入角
                        double phi = std::atan2(windSpeed * (1.0 - a[idx]), omega * r * (1.0 + ap[idx]));
                        // 计算迎角
                        double aoa = (phi - twist) * 180.0 / M_PI;
                        if (std::isnan(aoa) || std::abs(aoa) == 180.0)
                        {
                            aoa = 0.0;
                        }
                        perf.setAoaAt(b, t, i) = aoa;

                        int airfoilIdx = geom.airfoilIndex[i];
                        perf.setClAt(b, t, i) = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cl, aoa);
                        perf.setCdAt(b, t, i) = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cd, aoa);

                        double cl = perf.clAt(b, t, i);
                        double cd = perf.cdAt(b, t, i);

                        // 损失因子
                        double ftip = 2.0 / pi * std::acos(std::exp(-turbineParams.nBlades * (turbineParams.rTip - rNodes[i]) / (2.0 * rNodes[i] * std::sin(phi))));
                        double fhub = 2.0 / pi * std::acos(std::exp(-turbineParams.nBlades * (rNodes[i] - turbineParams.rHub) / (2.0 * turbineParams.rHub * std::sin(phi))));
                        double f = ftip * fhub;

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
                std::cout << "BEM converged after " << iter + 1 << " iterations" << std::endl;
                break;
            }
            if (iter == maxIterBEM - 1)
            {
                std::cerr << "Warning: BEM failed to converge after " << maxIterBEM << " iterations" << std::endl;
            }
        }
    }

} // namespace fvw