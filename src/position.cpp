#include "position.h"
#include <cmath>

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

    Vec3 &PositionData::leadAt(int b, int t, int i)
    {
        return lead[b * nTimesteps * nTrail + t * nTrail + i];
    }

    Vec3 &PositionData::quarterAt(int b, int t, int i)
    {
        return quarter[b * nTimesteps * nTrail + t * nTrail + i];
    }

    Vec3 &PositionData::trailAt(int b, int t, int i)
    {
        return trail[b * nTimesteps * nTrail + t * nTrail + i];
    }

    Vec3 &PositionData::collocAt(int b, int t, int i)
    {
        return colloc[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PositionData::boundAt(int b, int t, int i)
    {
        return bound[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PositionData::endAt(int b, int t, int i)
    {
        return end[b * nTimesteps * nShed + t * nShed + i];
    }

    Vec3 &PositionData::hubAt(int t) { return hub[t]; }
    Vec3 &PositionData::platformAt(int t) { return platform[t]; }

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
            pos.hubAt(t) = hubPos;
            pos.platformAt(t) = Vec3(0.0, 0.0, 0.0);
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

        // Rotation sequence
        std::string rseq = "zzzyxxyzxyz";
        std::vector<double> hubRotationSequence(4, 0.0);
        std::vector<double> bladePitch(simParams.timesteps, 0.0);

        // Compute positions
        for (int b = 0; b < turbineParams.nBlades; ++b)
        {
            std::vector<double> bladeRotseq = {90.0, bladePitch[0], 0.0, aziOri[b], 0.0, 0.0,
                                               hubRotationSequence[0], hubRotationSequence[1],
                                               hubRotationSequence[2], hubRotationSequence[3]};

            for (int t = 0; t < simParams.timesteps; ++t)
            {
                bladeRotseq[4] = azimuth[t]; // Update azimuth
                for (int i = 0; i < geom.lead.size(); ++i)
                {
                    std::vector<double> totalRotseq = {geom.twistTrailing[i]};
                    totalRotseq.insert(totalRotseq.end(), bladeRotseq.begin(), bladeRotseq.end());
                    pos.leadAt(b, t, i) = DCMRot(geom.lead[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.quarterAt(b, t, i) = DCMRot(geom.quarter[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.trailAt(b, t, i) = DCMRot(geom.trail[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                }
                for (int i = 0; i < geom.colloc.size(); ++i)
                {
                    std::vector<double> totalRotseq = {geom.twistShedding[i]};
                    totalRotseq.insert(totalRotseq.end(), bladeRotseq.begin(), bladeRotseq.end());
                    pos.collocAt(b, t, i) = DCMRot(geom.colloc[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.boundAt(b, t, i) = DCMRot(geom.bound[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                    pos.endAt(b, t, i) = DCMRot(geom.end[i], totalRotseq, {}, rseq, 0) + pos.hubAt(t);
                }
            }
        }
    }

} // namespace fvw