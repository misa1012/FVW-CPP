#ifndef FVW_CONFIG_H
#define FVW_CONFIG_H

#include <string>
#include <vector>
#include "core/geometry.h"

namespace fvw {

    struct PerturbationConfig {
        std::string name;
        PerturbationType type;
        double amplitude_deg;
        double freqFactor; // For dynamic cases: freq = freqFactor * rotationalFreq
    };

    struct GlobalConfig {
        TurbineParams turbine;
        SimParams sim; // Base simulation parameters
        std::vector<PerturbationConfig> perturbations;
        std::string caseName;   // Optional override for output case name
        std::string outputRoot; // Optional override for results root
    };

    class ConfigLoader {
    public:
        // Load configuration from a JSON file
        static GlobalConfig load(const std::string& filepath);
    };

} // namespace fvw

#endif // FVW_CONFIG_H
