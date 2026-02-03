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

    // Read raw config text for traceability
    std::string config_text;
    {
        std::ifstream cfg_in(config_path);
        if (cfg_in) {
            std::ostringstream ss;
            ss << cfg_in.rdbuf();
            config_text = ss.str();
        }
    }

    // 3. Setup Global Output
    std::filesystem::path output_path;
    std::filesystem::path exe_path = std::filesystem::canonical("/proc/self/exe");
    if (argc > 2) {
        output_path = argv[2];
    } else if (!config.outputRoot.empty()) {
        output_path = config.outputRoot;
    } else {
        std::filesystem::path project_root = exe_path.parent_path().parent_path(); // build/fvw_cpp -> FVW-CPP
        output_path = project_root / "results";
    }
    std::filesystem::create_directories(output_path);
    const std::string rootOutput = output_path.string();

    // 4. Run Cases
    auto batch_start = std::chrono::high_resolution_clock::now();
    
    std::vector<std::string> completed_cases;

    // Build run metadata (shared across all cases)
    fvw::RunMetadata runMeta;
    runMeta.config_path = config_path;
    runMeta.config_text = config_text;
    runMeta.exe_path = exe_path.string();
    {
        std::ostringstream args;
        for (int i = 0; i < argc; ++i) {
            if (i > 0) args << " ";
            args << argv[i];
        }
        runMeta.cli_args = args.str();
    }

    if (config.perturbations.empty()) {
        fvw::cli::print_warning("No perturbations defined. Running default baseline.");
        fvw::PerturbationConfig default_pc;
        default_pc.name = "default_baseline";
        default_pc.type = fvw::PerturbationType::None;
        
        if (!config.caseName.empty()) {
            default_pc.name = config.caseName;
        }
        fvw::SimulationRunner runner(default_pc, config, rootOutput, runMeta);
        runner.initialize();
        runner.run();
        runner.finalize();
        completed_cases.push_back("default_baseline");
    } else {
        std::cout << fvw::cli::BOLD << "Starting batch of " << config.perturbations.size() << " cases." << fvw::cli::RESET << "\n";
        for (const auto &pc : config.perturbations)
        {
            fvw::PerturbationConfig pc_case = pc;
            if (!config.caseName.empty()) {
                if (config.perturbations.size() == 1) {
                    pc_case.name = config.caseName;
                } else {
                    pc_case.name = config.caseName + "_" + pc.name;
                }
            }
            fvw::SimulationRunner runner(pc_case, config, rootOutput, runMeta);
            runner.initialize();
            runner.run();
            runner.finalize();
            completed_cases.push_back(pc_case.name);
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
