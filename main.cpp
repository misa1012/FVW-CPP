#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>

#include "simulation/simulation_runner.h"
#include "io/cli_utils.h"

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
    const std::string rootOutput = "../results";
    std::filesystem::create_directories(rootOutput);

    // 4. Run Cases
    const bool projectToGrid = false; 
    const bool computeProbes = false; 

    auto batch_start = std::chrono::high_resolution_clock::now();
    
    std::vector<std::string> completed_cases;

    if (config.perturbations.empty()) {
        fvw::cli::print_warning("No perturbations defined. Running default baseline.");
        fvw::PerturbationConfig default_pc;
        default_pc.name = "default_baseline";
        default_pc.type = fvw::PerturbationType::None;
        
        fvw::SimulationRunner runner(default_pc, config, rootOutput);
        runner.initialize();
        runner.run();
        runner.finalize(projectToGrid, computeProbes);
        completed_cases.push_back("default_baseline");
    } else {
        std::cout << fvw::cli::BOLD << "Starting batch of " << config.perturbations.size() << " cases." << fvw::cli::RESET << "\n";
        for (const auto &pc : config.perturbations)
        {
            fvw::SimulationRunner runner(pc, config, rootOutput);
            runner.initialize();
            runner.run();
            runner.finalize(projectToGrid, computeProbes);
            completed_cases.push_back(pc.name);
        }
    }
    
    auto batch_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::milliseconds>(batch_end - batch_start).count() / 1000.0;

    fvw::cli::print_header("Batch Summary");
    std::cout << fvw::cli::BOLD << std::left << std::setw(30) << "Case Name" << "Status" << fvw::cli::RESET << "\n";
    std::cout << "------------------------------------------\n";
    for(const auto& name : completed_cases) {
        std::cout << std::left << std::setw(30) << name << fvw::cli::GREEN << "COMPLETED" << fvw::cli::RESET << "\n";
    }
    std::cout << "\n" << fvw::cli::BOLD << "Total Time: " << total_time << " s" << fvw::cli::RESET << "\n";

    return 0;
}
