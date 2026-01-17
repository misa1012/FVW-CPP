#include "core/position.h"
#include <cmath>
#include <iostream>

namespace fvw
{

    PositionData::PositionData(int nBlades_, int timesteps_, int nTrail_, int nShed_)
        : nBlades(nBlades_), nTimesteps(timesteps_), nTrail(nTrail_), nShed(nShed_)
    {
        lead.resize(nBlades * nTimesteps * nTrail);
        quarter.resize(nBlades * nTimesteps * nTrail);
        trail.resize(nBlades * nTimesteps * nTrail);
        colloc.resize(nBlades * nTimesteps * nShed);
        bound.resize(nBlades * nTimesteps * nShed);
        end.resize(nBlades * nTimesteps * nShed);
        hub.resize(nTimesteps);
        platform.resize(nTimesteps);
    }

    // 只读访问器
    const Vec3 &PositionData::leadAt(int b, int t, int i) const
    {
        return lead[b * nTimesteps * nTrail + t * nTrail + i];
    }

    const Vec3 &PositionData::quarterAt(int b, int t, int i) const
    {
        return quarter[b * nTimesteps * nTrail + t * nTrail + i];
    }

    const Vec3 &PositionData::trailAt(int b, int t, int i) const
    {
        return trail[b * nTimesteps * nTrail + t * nTrail + i];
    }

    const Vec3 &PositionData::collocAt(int b, int t, int i) const
    {
        return colloc[b * nTimesteps * nShed + t * nShed + i];
    }

    const Vec3 &PositionData::boundAt(int b, int t, int i) const
    {
        return bound[b * nTimesteps * nShed + t * nShed + i];
    }

    const Vec3 &PositionData::endAt(int b, int t, int i) const
    {
        return end[b * nTimesteps * nShed + t * nShed + i];
    }

    const Vec3 &PositionData::hubAt(int t) const
    {
        return hub[t];
    }

    const Vec3 &PositionData::platformAt(int t) const
    {
        return platform[t];
    }

    // 可写访问器
    Vec3 &PositionData::setLeadAt(int b, int t, int i)
    {
        return lead[b * nTimesteps * nTrail + t * nTrail + i];
    }

    Vec3 &PositionData::setQuarterAt(int b, int t, int i)
    {
        return quarter[b * nTimesteps * nTrail + t * nTrail + i];
    }

    Vec3 &PositionData::setTrailAt(int b, int t, int i)
    {
        return trail[b * nTimesteps * nTrail + t * nTrail + i];
    }

    Vec3 &PositionData::setCollocAt(int b, int t, int i)
    {
        return colloc[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PositionData::setBoundAt(int b, int t, int i)
    {
        return bound[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PositionData::setEndAt(int b, int t, int i)
    {
        return end[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PositionData::setHubAt(int t)
    {
        return hub[t];
    }

    Vec3 &PositionData::setPlatformAt(int t)
    {
        return platform[t];
    }

    Vec3 DCMRot(const Vec3 &x, const std::vector<double> &t,
                const std::vector<double> &A_init, const std::string &rotseq, int rev)
    {
        std::vector<double> A = A_init.empty() ? std::vector<double>{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0} : A_init;

        for (size_t c1 = 0; c1 < rotseq.size() && c1 < t.size(); ++c1)
        {
            double rad = t[c1] * M_PI / 180.0;
            double sint = std::sin(rad);
            double cost = std::cos(rad);
            std::vector<double> R(9);

            if (rotseq[c1] == 'x' || rotseq[c1] == 'X')
            {
                R = {1.0, 0.0, 0.0, 0.0, cost, -sint, 0.0, sint, cost};
            }
            else if (rotseq[c1] == 'y' || rotseq[c1] == 'Y')
            {
                R = {cost, 0.0, sint, 0.0, 1.0, 0.0, -sint, 0.0, cost};
            }
            else if (rotseq[c1] == 'z' || rotseq[c1] == 'Z')
            {
                R = {cost, -sint, 0.0, sint, cost, 0.0, 0.0, 0.0, 1.0};
            }

            std::vector<double> B(9);
            B[0] = R[0] * A[0] + R[1] * A[3] + R[2] * A[6];
            B[1] = R[0] * A[1] + R[1] * A[4] + R[2] * A[7];
            B[2] = R[0] * A[2] + R[1] * A[5] + R[2] * A[8];
            B[3] = R[3] * A[0] + R[4] * A[3] + R[5] * A[6];
            B[4] = R[3] * A[1] + R[4] * A[4] + R[5] * A[7];
            B[5] = R[3] * A[2] + R[4] * A[5] + R[5] * A[8];
            B[6] = R[6] * A[0] + R[7] * A[3] + R[8] * A[6];
            B[7] = R[6] * A[1] + R[7] * A[4] + R[8] * A[7];
            B[8] = R[6] * A[2] + R[7] * A[5] + R[8] * A[8];
            A = B;
        }

        if (rev == 1)
        {
            std::vector<double> B = A;
            A[1] = B[3];
            A[2] = B[6];
            A[3] = B[1];
            A[5] = B[7];
            A[6] = B[2];
            A[7] = B[5];
        }

        Vec3 y;
        y.x = A[0] * x.x + A[1] * x.y + A[2] * x.z;
        y.y = A[3] * x.x + A[4] * x.y + A[5] * x.z;
        y.z = A[6] * x.x + A[7] * x.y + A[8] * x.z;
        return y;
    }

    void computePositions(PositionData &pos, const SimParams &simParams,
                          const TurbineParams &turbineParams, const BladeGeometry &geom)
    {
        // Hub and platform
        Vec3 hubPos(0.0, 0.0, 90.0);
        for (int t = 0; t < simParams.timesteps; ++t)
        {
            pos.setHubAt(t) = hubPos;
            pos.setPlatformAt(t) = Vec3(0.0, 0.0, 0.0);
        }

        // Azimuth
        std::vector<double> azimuth(simParams.timesteps);
        for (int t = 0; t < simParams.timesteps; ++t)
        {
            azimuth[t] = turbineParams.omega * (180.0 / M_PI) * (t * simParams.dt);
        }

        // Blade initial angles
        double aziStep = 360.0 / turbineParams.nBlades;
        std::vector<double> aziOri(turbineParams.nBlades);
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            aziOri[b] = b * aziStep;
        }

        // 创建一个基础的、不含扰动的桨距角数组 (比如全为0)
        std::vector<std::vector<double>> blade_pitches(simParams.timesteps, std::vector<double>(turbineParams.nBlades, 0.0));

        // 如果开启了扰动，则在基础值上添加正弦波动
        switch (simParams.perturbation.type)
        {
        case PerturbationType::CollectivePitch:
        {
            // 注意：这里是集体变桨，所以所有叶片在同一时刻的桨距角增量是相同的
            std::cout << "--- Applying Collective Pitch Perturbation ---" << std::endl;
            std::cout << "  - Amplitude: " << simParams.perturbation.amplitude_deg << " deg" << std::endl;
            std::cout << "  - Frequency: " << simParams.perturbation.frequency_hz << " Hz" << std::endl;

            double omega_perturb = 2.0 * M_PI * simParams.perturbation.frequency_hz;

            for (int t = 0; t < simParams.timesteps; ++t)
            {
                double currentTime = t * simParams.dt;
                double pitch_perturbation = simParams.perturbation.amplitude_deg * std::sin(omega_perturb * currentTime);

                // 将扰动量加到基础桨距角上 (这里假设基础值为0)
                // 其中，DCM的输入为度
                for (int b = 0; b < turbineParams.nBlades; ++b)
                {
                    blade_pitches[t][b] = pitch_perturbation;
                }
            }
        }
        break;
        case PerturbationType::AsymmetricStaticPitch:
        {
            std::cout << "--- Applying Asymmetric Static Pitch Perturbation ---" << std::endl;
            std::cout << "  - Delta Alpha: " << simParams.perturbation.amplitude_deg << " deg" << std::endl;

            std::vector<double> pitch_offsets(turbineParams.nBlades);

            // 假设3个叶片：一个为0，一个为+delta, 一个为-delta
            if (turbineParams.nBlades == 3)
            {
                pitch_offsets[0] = 0.0;
                pitch_offsets[1] = simParams.perturbation.amplitude_deg;
                pitch_offsets[2] = -simParams.perturbation.amplitude_deg;
            }
            else
            {
                std::cerr << "[Error] Blade number mismatch: Asymmetric Static Pitch requires 3 blades." << std::endl;
                for (int b = 0; b < turbineParams.nBlades; ++b)
                {
                    pitch_offsets[b] = simParams.perturbation.amplitude_deg * std::sin(b * 2.0 * M_PI / turbineParams.nBlades);
                }
            }

            // 将这个固定的偏移量应用到所有时间步
            for (int t = 0; t < simParams.timesteps; ++t)
            {
                for (int b = 0; b < turbineParams.nBlades; ++b)
                {
                    blade_pitches[t][b] = pitch_offsets[b];
                }
            }
            break;
        }
        case PerturbationType::None:
        default:
            // 默认情况下，所有桨距角都为0，无需操作
            break;
        }

        // Rotation sequence
        std::string rseq = "zzzyxxyzxyz";
        std::vector<double> hubRotationSequence(4, 0.0);

        // Compute positions
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            for (int t = 0; t < simParams.timesteps; ++t)
            {
                std::vector<double> bladeRotseq = {90.0, blade_pitches[t][b], 0.0, aziOri[b], azimuth[t], 0.0,
                                                   hubRotationSequence[0], hubRotationSequence[1],
                                                   hubRotationSequence[2], hubRotationSequence[3]};

                for (size_t i = 0; i < geom.lead.size(); ++i)
                {
                    std::vector<double> totalRotseq = {geom.twistTrailing[i]};
                    totalRotseq.insert(totalRotseq.end(), bladeRotseq.begin(), bladeRotseq.end());
                    pos.setLeadAt(b, t, i) = DCMRot(geom.lead[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.setQuarterAt(b, t, i) = DCMRot(geom.quarter[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.setTrailAt(b, t, i) = DCMRot(geom.trail[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                }
                for (size_t i = 0; i < geom.colloc.size(); ++i)
                {
                    std::vector<double> totalRotseq = {geom.twistShedding[i]};
                    totalRotseq.insert(totalRotseq.end(), bladeRotseq.begin(), bladeRotseq.end());
                    pos.setCollocAt(b, t, i) = DCMRot(geom.colloc[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.setBoundAt(b, t, i) = DCMRot(geom.bound[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.setEndAt(b, t, i) = DCMRot(geom.end[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                }
            }
        }
    }

} // namespace fvw