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
        if (s == "VanGarrelUnitConsistent") return VortexCoreType::VanGarrelUnitConsistent;
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

    // Helper to map string to SegmentDistribution
    SegmentDistribution parseSegmentDistribution(const std::string& s) {
        if (s == "Linear") return SegmentDistribution::Linear;
        if (s == "Cosine") return SegmentDistribution::Cosine;
        throw std::runtime_error("Unknown SegmentDistribution: " + s);
    }
}

GlobalConfig ConfigLoader::load(const std::string& filepath) {
    GlobalConfig config;
    
    // Parse JSON file
    fvw::json::Value root = fvw::json::parse_file(filepath);
    bool log_verbose = false;
    if (root.contains("simulation") && root["simulation"].contains("logVerbose")) {
        log_verbose = root["simulation"]["logVerbose"].as_bool();
    }
    if (root.contains("caseName")) {
        config.caseName = root["caseName"].as_string();
    }
    if (root.contains("outputRoot")) {
        config.outputRoot = root["outputRoot"].as_string();
    }

    // --- Turbine Params ---
    if (!root.contains("turbine")) throw std::runtime_error("Missing 'turbine' section in config");
    const auto& turbineJson = root["turbine"];

    // 1. Identify Model (strict: must be provided)
    if (!turbineJson.contains("model")) {
        throw std::runtime_error("Missing required field: turbine.model");
    }
    std::string modelName = turbineJson["model"].as_string();
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
        if (log_verbose) {
            std::cout << "Loading default turbine parameters from: " << dataPath << std::endl;
        }
        try {
            fvw::json::Value defaults = fvw::json::parse_file(dataPath);
            
            // Populate defaults
            if (defaults.contains("rTip")) config.turbine.rTip = defaults["rTip"].as_double();
            if (defaults.contains("rHub")) config.turbine.rHub = defaults["rHub"].as_double();
            if (defaults.contains("hubHeight")) {
                config.turbine.hubHeight = defaults["hubHeight"].as_double();
                if (log_verbose) {
                    std::cout << "  - Loaded Default Hub Height: " << config.turbine.hubHeight << " m" << std::endl;
                }
            }
            if (defaults.contains("nBlades")) config.turbine.nBlades = defaults["nBlades"].as_int();

        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to load default params from " << dataPath << ": " << e.what() << std::endl;
        }
    }

    // 3. Overwrite with Config Values (Explicit Override)
    if (!turbineJson.contains("windSpeed")) throw std::runtime_error("Missing required field: turbine.windSpeed");
    if (!turbineJson.contains("rho")) throw std::runtime_error("Missing required field: turbine.rho");
    if (!turbineJson.contains("nSegments")) throw std::runtime_error("Missing required field: turbine.nSegments");
    if (!turbineJson.contains("segmentDistribution")) throw std::runtime_error("Missing required field: turbine.segmentDistribution");
    if (!turbineJson.contains("tsr")) throw std::runtime_error("Missing required field: turbine.tsr");

    if (turbineJson.contains("rHub")) config.turbine.rHub = turbineJson["rHub"].as_double();
    if (turbineJson.contains("rTip")) config.turbine.rTip = turbineJson["rTip"].as_double();
    if (turbineJson.contains("nBlades")) config.turbine.nBlades = turbineJson["nBlades"].as_int();

    config.turbine.windSpeed = turbineJson["windSpeed"].as_double();
    config.turbine.rho = turbineJson["rho"].as_double();
    config.turbine.nSegments = turbineJson["nSegments"].as_int();
    config.turbine.segmentDistribution = parseSegmentDistribution(turbineJson["segmentDistribution"].as_string());
    config.turbine.tsr = turbineJson["tsr"].as_double();
    
    // Validate critical parameters
    if (config.turbine.rTip <= 0.0) throw std::runtime_error("Critical parameter 'rTip' not set/invalid. Check turbine_params.json.");
    if (config.turbine.rHub <= 0.0) throw std::runtime_error("Critical parameter 'rHub' not set/invalid. Check turbine_params.json.");
    if (config.turbine.nBlades <= 0) throw std::runtime_error("Critical parameter 'nBlades' not set/invalid.");
    if (config.turbine.windSpeed <= 0.0) throw std::runtime_error("Critical parameter 'windSpeed' not set/invalid.");
    if (config.turbine.rho <= 0.0) throw std::runtime_error("Critical parameter 'rho' not set/invalid.");
    if (config.turbine.nSegments <= 0) throw std::runtime_error("Critical parameter 'nSegments' not set/invalid.");
    if (config.turbine.tsr <= 0.0) throw std::runtime_error("Critical parameter 'tsr' not set/invalid.");
    
    // Strict requirement: hubHeight must come from data file (or default loaded from there)
    // We initialized it to -1.0. If it wasn't loaded from defaults, it's invalid.
    if (config.turbine.hubHeight <= 0.0) {
        throw std::runtime_error("Critical parameter 'hubHeight' not set. It must be defined in turbine_params.json."); 
    }

    // Calculate derived parameter omega
    config.turbine.omega = config.turbine.tsr * config.turbine.windSpeed / config.turbine.rTip;

    // --- Simulation Params ---
    if (!root.contains("simulation")) throw std::runtime_error("Missing 'simulation' section in config");
    const auto& simJson = root["simulation"];
    
    // Check for revolution-based parameters first
    const bool has_steps = simJson.contains("stepsPerRevolution");
    const bool has_revs = simJson.contains("numRevolutions");
    const bool has_dt = simJson.contains("dt");
    const bool has_total = simJson.contains("totalTime");

    if (has_steps && has_revs) {
        config.sim.stepsPerRevolution = simJson["stepsPerRevolution"].as_int();
        config.sim.numRevolutions = simJson["numRevolutions"].as_double();
        
        // Calculate period of one revolution
        double period = 2.0 * M_PI / config.turbine.omega;
        
        // Calculate dt and totalTime from revolution-based parameters
        double raw_dt = period / config.sim.stepsPerRevolution;
        // Enforce fixed precision (6 decimal places) as requested
        config.sim.dt = std::round(raw_dt * 1000000.0) / 1000000.0;
        config.sim.totalTime = config.sim.numRevolutions * period;
        
        if (log_verbose) {
            std::cout << "Using revolution-based parameters: " 
                      << config.sim.stepsPerRevolution << " steps/rev, "
                      << config.sim.numRevolutions << " revolutions" << std::endl;
            std::cout << "Calculated: dt = " << config.sim.dt << " s, totalTime = " 
                      << config.sim.totalTime << " s" << std::endl;
        }
    } else {
        // Fallback to direct specification
        if (!(has_dt && has_total)) {
            throw std::runtime_error("Missing required simulation timing fields: either (stepsPerRevolution & numRevolutions) or (dt & totalTime).");
        }
        config.sim.dt = simJson["dt"].as_double();
        config.sim.totalTime = simJson["totalTime"].as_double();
    }
    if (!simJson.contains("outputFrequency")) throw std::runtime_error("Missing required field: simulation.outputFrequency");
    if (!simJson.contains("cutoffParam")) throw std::runtime_error("Missing required field: simulation.cutoffParam");
    if (!simJson.contains("coreType")) throw std::runtime_error("Missing required field: simulation.coreType");
    if (!simJson.contains("vortexModel")) throw std::runtime_error("Missing required field: simulation.vortexModel");

    config.sim.outputFrequency = simJson["outputFrequency"].as_int();
    if (simJson.contains("probeFrequency")) config.sim.probeFrequency = simJson["probeFrequency"].as_int();
    if (simJson.contains("computeProbes")) config.sim.computeProbes = simJson["computeProbes"].as_bool();
    config.sim.cutoffParam = simJson["cutoffParam"].as_double();
    config.sim.coreType = parseCoreType(simJson["coreType"].as_string());
    config.sim.vortexModel = parseVortexModel(simJson["vortexModel"].as_string());

    // BEM Solver Settings (Strict)
    if (!simJson.contains("bemTolerance")) throw std::runtime_error("Missing required field: simulation.bemTolerance");
    if (!simJson.contains("bemMaxIterations")) throw std::runtime_error("Missing required field: simulation.bemMaxIterations");
    if (!simJson.contains("bemRelaxation")) throw std::runtime_error("Missing required field: simulation.bemRelaxation");
    config.sim.bemTolerance = simJson["bemTolerance"].as_double();
    config.sim.bemMaxIterations = simJson["bemMaxIterations"].as_int();
    config.sim.bemRelaxation = simJson["bemRelaxation"].as_double();

    // Kutta-Joukowski Iteration Settings (Optional)
    if (simJson.contains("kuttaTolerance")) config.sim.kuttaTolerance = simJson["kuttaTolerance"].as_double();
    if (simJson.contains("kuttaMaxIterations")) config.sim.kuttaMaxIterations = simJson["kuttaMaxIterations"].as_int();
    if (simJson.contains("kuttaRelaxation")) config.sim.kuttaRelaxation = simJson["kuttaRelaxation"].as_double();

    // Optional logging settings
    if (simJson.contains("logStepTiming")) config.sim.logStepTiming = simJson["logStepTiming"].as_bool();
    if (simJson.contains("logVerbose")) config.sim.logVerbose = simJson["logVerbose"].as_bool();
    if (simJson.contains("logPerf")) config.sim.logPerf = simJson["logPerf"].as_bool();

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
