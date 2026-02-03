#include "core/geometry.h"
#include "core/utils.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace fvw
{
    BladeDefinition loadBladeDefinition(const std::string &csvPath) {
        BladeDefinition bd;
        std::ifstream file(csvPath);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open blade geometry file: " << csvPath << std::endl;
            // Return empty or throw, here we return empty which will likely cause issues downstream but handled simply for now
            return bd;
        }

        std::string line;
        // Skip header
        std::getline(file, line);

        while (std::getline(file, line)) {
            if (line.empty()) continue;
            // Handle simple CSV parsing
            std::replace(line.begin(), line.end(), ',', ' '); 
            std::istringstream iss(line);
            double r, c, t;
            int idx;
            if (iss >> r >> c >> t >> idx) {
                bd.r.push_back(r);
                bd.chord.push_back(c);
                bd.twist.push_back(t);
                bd.airfoilIndex.push_back(idx);
            }
        }
        return bd;
    }

    BladeGeometry computeBladeGeometry(const TurbineParams &params, const BladeDefinition &rawDist)
    {
        BladeGeometry geom;
        int nTrail = params.nSegments + 1;
        int nShed = params.nSegments;

        const auto& originalR = rawDist.r;
        const auto& originalChord = rawDist.chord;
        const auto& originalTwist = rawDist.twist;
        const auto& originalNFoil = rawDist.airfoilIndex;

        if (originalR.empty()) {
             std::cerr << "Error: Blade definition is empty!" << std::endl;
             return geom;
        }

        // Cosine distribution for trailing nodes
        geom.rTrailing.resize(nTrail);
        geom.chordTrailing.resize(nTrail);
        geom.twistTrailing.resize(nTrail);
        std::vector<double> theta(nTrail);
        for (int i = 0; i < nTrail; ++i)
        {
            theta[i] = M_PI * i / (nTrail - 1);
            geom.rTrailing[i] = params.rHub + (params.rTip - params.rHub) * (1 - std::cos(theta[i])) / 2;
            geom.chordTrailing[i] = interpolate(originalR, originalChord, geom.rTrailing[i]);
            geom.twistTrailing[i] = interpolate(originalR, originalTwist, geom.rTrailing[i]);
        }

        // Shedding nodes
        geom.rShedding.resize(nShed);
        geom.chordShedding.resize(nShed);
        geom.twistShedding.resize(nShed);
        geom.airfoilIndex.resize(nShed);
        for (int i = 0; i < nShed; ++i)
        {
            geom.rShedding[i] = (geom.rTrailing[i] + geom.rTrailing[i + 1]) / 2;
            geom.chordShedding[i] = interpolate(geom.rTrailing, geom.chordTrailing, geom.rShedding[i]);
            geom.twistShedding[i] = interpolate(geom.rTrailing, geom.twistTrailing, geom.rShedding[i]);
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