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

        // Reading physical parameters
        double windSpeed = turbineParams.windSpeed;
        double omega = turbineParams.omega;
        const double pi = M_PI;

        // --- 1D Solver storage (for one blade, one timestep) ---
        int nShed = perf.getShed();
        std::vector<double> a(nShed, 0.0);
        std::vector<double> ap(nShed, 0.0);
        std::vector<double> a0(nShed);
        std::vector<double> ap0(nShed);
        std::vector<double> solidity(nShed);

        // Pre-calculate solidity
        for (int i = 0; i < nShed; ++i)
        {
            solidity[i] = (turbineParams.nBlades * geom.chordShedding[i]) / (2 * pi * geom.rShedding[i]);
        }

        // --- Initial Guess (Vectorized for segments) ---
        for (int i = 0; i < nShed; ++i)
        {
            double lambda_r = omega * geom.rShedding[i] / windSpeed;
            double twist = -geom.twistShedding[i] * pi / 180.0;
            double expr = 4.0 - 4.0 * pi * lambda_r * solidity[i] +
                          pi * lambda_r * lambda_r * solidity[i] * (8.0 * twist + pi * solidity[i]);
            
            if (expr < 0.0 || std::isnan(expr))
            {
                a[i] = 0.33;
                // Only log warning once per segment if needed, or suppress
                 std::cerr << "[Warning] Invalid sqrt_expr at i=" << i << ", expr=" << expr << ", set a=0.33" << std::endl;
            }
            else
            {
                a[i] = 0.25 * (2.0 + pi * lambda_r * solidity[i] - std::sqrt(expr));
            }
            ap[i] = 0.0;
            
            if (std::isnan(a[i]) || std::isinf(a[i]))
            {
                a[i] = 0.33;
                std::cerr << "[Warning] Initial a NaN/Inf at i=" << i << ", set a=0.33" << std::endl;
            }
        }

        // --- BEM Iteration (Vectorized for segments) ---
        // We iterate for all segments together until ALL converge (simplest approach for vectorization)
        // Or strictly, we keep iterating.
        
        // Storage for results to broadcast later
        std::vector<double> final_aoa(nShed);
        std::vector<double> final_cl(nShed);
        std::vector<double> final_cd(nShed);

        for (int iter = 0; iter < maxIterBEM; ++iter)
        {
            bool global_converged = true;

            for (int i = 0; i < nShed; ++i)
            {
                // Save old values
                a0[i] = a[i];
                ap0[i] = ap[i];

                double r = geom.rShedding[i];
                double twist = -1 * geom.twistShedding[i] * M_PI / 180.0; 

                // Sanity check
                if (std::isnan(a[i]) || std::isnan(ap[i])) {
                    a[i] = 0.33; ap[i] = 0.0;
                }

                double phi = std::atan2(windSpeed * (1.0 - a[i]), omega * r * (1.0 + ap[i]));
                double aoa = (phi - twist) * 180.0 / M_PI;

                if (std::isnan(aoa) || std::abs(aoa) == 180.0) aoa = 0.0;

                // Lookup airfoil
                int airfoilIdx = geom.airfoilIndex[i];
                double cl = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cl, aoa);
                double cd = interpolate(airfoils[airfoilIdx].aoa, airfoils[airfoilIdx].cd, aoa);
                
                // Store temporarily
                final_aoa[i] = aoa;
                final_cl[i] = cl;
                final_cd[i] = cd;

                // Tip/Hub Loss
                double ftip = 2.0 / pi * std::acos(std::exp(-turbineParams.nBlades * (turbineParams.rTip - geom.rShedding[i]) / (2.0 * geom.rShedding[i] * std::sin(phi))));
                double fhub = 2.0 / pi * std::acos(std::exp(-turbineParams.nBlades * (geom.rShedding[i] - turbineParams.rHub) / (2.0 * turbineParams.rHub * std::sin(phi))));
                double f = ftip * fhub;
                if (std::isnan(f) || f <= 0.0) f = 1.0; 

                // Thrust Coeff
                double ct = solidity[i] * (1.0 - a[i]) * (1.0 - a[i]) *
                            (cl * std::cos(phi) + cd * std::sin(phi)) /
                            (std::sin(phi) * std::sin(phi));

                // Update a
                double denom = solidity[i] * (cl * std::cos(phi) + cd * std::sin(phi));
                double a_new = denom > 1e-10 ? 1.0 / (1.0 + 4.0 * f * std::sin(phi) * std::sin(phi) / denom) : a[i];
                
                // Glauert correction
                if (ct > 0.96 * f)
                {
                    a_new = (18.0 * f - 20.0 - 3.0 * std::sqrt(ct * (50.0 - 36.0 * f) + 12.0 * f * (3.0 * f - 4.0))) /
                            (36.0 * f - 50.0);
                }

                // Update ap
                double ap_new = denom > 1e-10 ? 1.0 / (4.0 * f * std::cos(phi) * std::sin(phi) /
                                                           (solidity[i] * (cl * std::sin(phi) - cd * std::cos(phi))) -
                                                       1.0)
                                              : ap[i];

                // Relaxation
                double da = std::abs(a_new - a0[i]);
                double dap = std::abs(ap_new - ap0[i]);

                a[i] = (1.0 - weightFactor) * a[i] + weightFactor * a_new;
                ap[i] = (1.0 - weightFactor) * ap[i] + weightFactor * ap_new;

                if (da > tolBEM || dap > tolBEM) global_converged = false;
            }

            if (global_converged) {
                if (simParams.logVerbose) {
                    Logger::log("SOLVER", "BEM Converged: " + std::to_string(iter + 1) + " iters");
                }
                break;
            }
            if (iter == maxIterBEM - 1) {
                std::cerr << "Warning: BEM failed to converge after " << maxIterBEM << " iterations" << std::endl;
            }
        }


        // --- Broadcast Results to PerformanceData (All blades, All timesteps) ---
        for (int b = 0; b < perf.getBlades(); ++b)
        {
            for (int t = 0; t < perf.getTimesteps(); ++t)
            {
                for (int i = 0; i < nShed; ++i)
                {
                    perf.setAoaAt(b, t, i) = final_aoa[i];
                    perf.setClAt(b, t, i) = final_cl[i];
                    perf.setCdAt(b, t, i) = final_cd[i];
                    
                    // We can also store induced velocities if PerformanceData requires them, 
                    // but standard init focuses on AoA/Cl/Cd for the first wake step.
                    // If bound gamma is needed:
                    // Gamma = 0.5 * W * c * Cl (This logic might be inside KuttaJoukowski, 
                    // but if needed here we can compute it if we knew W, relative velocity).
                    // Typically BEM init sets the "Target" AoA/Cl/Cd.
                }
            }
        }
    }

} // namespace fvw
