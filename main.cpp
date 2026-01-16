#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>

#include "simulation_runner.h"

// Clean output dir (moved to SimulationRunner but good to have helper or use from class?)
// Actually SimulationRunner handles its own output directory.
// We just need main to drive the batch.

int main(int argc, char *argv[]) {
    // 1. Determine config file path
    std::string config_path = "config.json";
    if (argc > 1) {
        config_path = argv[1];
    }
    
    // Check if configuration exists
    if (!std::filesystem::exists(config_path)) {
        if (std::filesystem::exists("../" + config_path)) {
             config_path = "../" + config_path;
        } else {
             std::cerr << "Error: Config file not found: " << config_path << "\n";
             return 1;
        }
    }
    
    std::cout << "Loading configuration from: " << config_path << std::endl;

    // 2. Load Configuration
    fvw::GlobalConfig config;
    try {
        config = fvw::ConfigLoader::load(config_path);
    } catch (const std::exception& e) {
        std::cerr << "Error loading configuration: " << e.what() << "\n";
        return 1;
    }

#ifdef NDEBUG
    std::cout << "This is a Release build." << std::endl;
#else
    std::cout << "This is a Debug build." << std::endl;
#endif

    // 3. Setup Global Output
    const std::string rootOutput = "./results";
    std::filesystem::create_directories(rootOutput);

    // 4. Run Cases
    // Control flags (could be exposed via CLI args later)
    const bool projectToGrid = false; 
    const bool computeProbes = false; 

    auto total_start = std::chrono::high_resolution_clock::now();
    
    if (config.perturbations.empty()) {
        std::cout << "No perturbations defined in config. Running single default baseline.\n";
        fvw::PerturbationConfig default_pc;
        default_pc.name = "default_baseline";
        default_pc.type = fvw::PerturbationType::None;
        
        fvw::SimulationRunner runner(default_pc, config, rootOutput);
        runner.initialize();
        runner.run();
        runner.finalize(projectToGrid, computeProbes);
    } else {
        std::cout << "Starting batch of " << config.perturbations.size() << " cases.\n";
        for (const auto &pc : config.perturbations)
        {
            fvw::SimulationRunner runner(pc, config, rootOutput);
            runner.initialize();
            runner.run();
            runner.finalize(projectToGrid, computeProbes);
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(total_end - total_start);

    std::cout << "\n[ALL DONE] Batch completed in "
              << total_duration.count() / 1e6 << " s" << std::endl;

    return 0;
}
