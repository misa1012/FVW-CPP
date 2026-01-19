#include "io/config.h"
#include "io/json_utils.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <filesystem>

namespace fvw {

namespace {
    // Helper to map string to VortexModelType
    VortexModelType parseVortexModel(const std::string& s) {
        if (s == "Constant") return VortexModelType::Constant;
        if (s == "GammaDecay") return VortexModelType::GammaDecay;
        throw std::runtime_error("Unknown VortexModelType: " + s);
    }

    // Helper to map string to VortexCoreType
    VortexCoreType parseCoreType(const std::string& s) {
        if (s == "VanGarrel") return VortexCoreType::VanGarrel;
        if (s == "ChordBasedCore") return VortexCoreType::ChordBasedCore;
        throw std::runtime_error("Unknown VortexCoreType: " + s);
    }

    // Helper to map string to PerturbationType
    PerturbationType parsePerturbationType(const std::string& s) {
        if (s == "None") return PerturbationType::None;
        if (s == "CollectivePitch") return PerturbationType::CollectivePitch;
        if (s == "AsymmetricStaticPitch") return PerturbationType::AsymmetricStaticPitch;
        throw std::runtime_error("Unknown PerturbationType: " + s);
    }
}

GlobalConfig ConfigLoader::load(const std::string& filepath) {
    GlobalConfig config;
    
    // Parse JSON file
    fvw::json::Value root = fvw::json::parse_file(filepath);

    // --- Turbine Params ---
    if (!root.contains("turbine")) throw std::runtime_error("Missing 'turbine' section in config");
    const auto& turbineJson = root["turbine"];

    // 1. Identify Model
    std::string modelName = "NREL_5MW"; // default fallback or required? 
    if (turbineJson.contains("model")) {
        modelName = turbineJson["model"].as_string();
    }
    config.turbine.model = modelName;

    // 2. Load Defaults from Data File if exists
    // Path resolution logic similar to simulation_runner
    std::string dataPath = "data/" + modelName + "/turbine_params.json";
    if (!std::filesystem::exists(dataPath)) {
        if (std::filesystem::exists("../" + dataPath)) {
            dataPath = "../" + dataPath;
        }
    }

    if (std::filesystem::exists(dataPath)) {
        std::cout << "Loading default turbine parameters from: " << dataPath << std::endl;
        try {
            fvw::json::Value defaults = fvw::json::parse_file(dataPath);
            
            // Populate defaults
            if (defaults.contains("rTip")) config.turbine.rTip = defaults["rTip"].as_double();
            if (defaults.contains("rHub")) config.turbine.rHub = defaults["rHub"].as_double();
            if (defaults.contains("hubHeight")) {
                config.turbine.hubHeight = defaults["hubHeight"].as_double();
                std::cout << "  - Loaded Default Hub Height: " << config.turbine.hubHeight << " m" << std::endl;
            }
            if (defaults.contains("nBlades")) config.turbine.nBlades = defaults["nBlades"].as_int();

        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to load default params from " << dataPath << ": " << e.what() << std::endl;
        }
    }

    // 3. Overwrite with Config Values (Explicit Override)
    if (turbineJson.contains("windSpeed")) config.turbine.windSpeed = turbineJson["windSpeed"].as_double();
    if (turbineJson.contains("rho")) config.turbine.rho = turbineJson["rho"].as_double();
    if (turbineJson.contains("rHub")) config.turbine.rHub = turbineJson["rHub"].as_double();
    if (turbineJson.contains("rTip")) config.turbine.rTip = turbineJson["rTip"].as_double();
    if (turbineJson.contains("hubHeight")) {
        config.turbine.hubHeight = turbineJson["hubHeight"].as_double();
        std::cout << "  - Configuration Overriding Hub Height: " << config.turbine.hubHeight << " m" << std::endl;
    }
    if (turbineJson.contains("nBlades")) config.turbine.nBlades = turbineJson["nBlades"].as_int();
    if (turbineJson.contains("nSegments")) config.turbine.nSegments = turbineJson["nSegments"].as_int();
    if (turbineJson.contains("tsr")) config.turbine.tsr = turbineJson["tsr"].as_double();
    
    // Validate required fields (if defaults didn't cover them)
    // windSpeed, nSegments, tsr are unlikely to be in defaults, so usually required.
    // If they were not set (initially 0 or garbage), we might want to check here.
    // However, C++ structs without constructor might have garbage.
    // Assuming user provides valid config.

    // Calculate derived parameter omega
    config.turbine.omega = config.turbine.tsr * config.turbine.windSpeed / config.turbine.rTip;

    // --- Simulation Params ---
    if (!root.contains("simulation")) throw std::runtime_error("Missing 'simulation' section in config");
    const auto& simJson = root["simulation"];
    
    // Check for revolution-based parameters first
    if (simJson.contains("stepsPerRevolution") && simJson.contains("numRevolutions")) {
        config.sim.stepsPerRevolution = simJson["stepsPerRevolution"].as_int();
        config.sim.numRevolutions = simJson["numRevolutions"].as_double();
        
        // Calculate period of one revolution
        double period = 2.0 * M_PI / config.turbine.omega;
        
        // Calculate dt and totalTime from revolution-based parameters
        config.sim.dt = period / config.sim.stepsPerRevolution;
        config.sim.totalTime = config.sim.numRevolutions * period;
        
        std::cout << "Using revolution-based parameters: " 
                  << config.sim.stepsPerRevolution << " steps/rev, "
                  << config.sim.numRevolutions << " revolutions" << std::endl;
        std::cout << "Calculated: dt = " << config.sim.dt << " s, totalTime = " 
                  << config.sim.totalTime << " s" << std::endl;
    } else {
        // Fallback to direct specification
        config.sim.dt = simJson["dt"].as_double();
        config.sim.totalTime = simJson["totalTime"].as_double();
    }
    config.sim.outputFrequency = simJson["outputFrequency"].as_int();
    config.sim.cutoffParam = simJson["cutoffParam"].as_double();
    config.sim.coreType = parseCoreType(simJson["coreType"].as_string());
    config.sim.vortexModel = parseVortexModel(simJson["vortexModel"].as_string());

    // BEM Solver Settings (Optional, with defaults)
    if (simJson.contains("bemTolerance")) config.sim.bemTolerance = simJson["bemTolerance"].as_double();
    if (simJson.contains("bemMaxIterations")) config.sim.bemMaxIterations = simJson["bemMaxIterations"].as_int();
    if (simJson.contains("bemRelaxation")) config.sim.bemRelaxation = simJson["bemRelaxation"].as_double();

    // Probe Settings (Optional)
    if (simJson.contains("probeFrequency")) config.sim.probeFrequency = simJson["probeFrequency"].as_int();
    else config.sim.probeFrequency = config.sim.outputFrequency; // Default to output frequency
    
    if (simJson.contains("computeProbes")) config.sim.computeProbes = simJson["computeProbes"].as_bool();
    else config.sim.computeProbes = false; // Default to false
    
    // Calculate derived timesteps (including t=0)
    config.sim.timesteps = static_cast<int>(config.sim.totalTime / config.sim.dt) + 1;

    // Default base perturbation to None
    config.sim.perturbation.type = PerturbationType::None;
    config.sim.perturbation.amplitude_deg = 0.0;
    config.sim.perturbation.frequency_hz = 0.0;

    // --- Perturbations List ---
    if (root.contains("perturbations")) {
        const auto& pertArray = root["perturbations"].as_array();
        for (const auto& item : pertArray) {
            PerturbationConfig pc;
            pc.name = item["name"].as_string();
            pc.type = parsePerturbationType(item["type"].as_string());
            pc.amplitude_deg = item["amplitude"].as_double(); // JSON key "amplitude" maps to amplitude_deg
            pc.freqFactor = item["freqFactor"].as_double();
            config.perturbations.push_back(pc);
        }
    }

    return config;
}

} // namespace fvw
