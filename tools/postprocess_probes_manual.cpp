#include <iostream>
#include <string>
#include <vector>
#include "io/postprocess.h"
#include "io/config.h"
#include "simulation/simulation_runner.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <config_file> [override_h5_file]\n";
        return 1;
    }

    std::string config_path = argv[1];
    
    // 1. Load Config
    std::cout << "Loading config from " << config_path << std::endl;
    fvw::GlobalConfig config;
    try {
        config = fvw::ConfigLoader::load(config_path);
    } catch (const std::exception& e) {
        std::cerr << "Error loading config: " << e.what() << std::endl;
        return 1;
    }

    // 2. Resolve Paths (Geometry)
    // We can reuse logic from SimulationRunner or manually replicate it
    // SimulationRunner runner(..., config, ...); runner.resolve_paths();
    // But instantiating runner might be heavy. Let's replicate simple logic or use a helper.
    // For now, let's just manually load geometry as we know it's NTNU.
    std::string geom_path = "data/" + config.turbine.model + "/blade_geometry.csv";
    std::cout << "Loading geometry from " << geom_path << std::endl;
    
    fvw::BladeDefinition bladeDef = fvw::loadBladeDefinition(geom_path);
    if (bladeDef.r.empty()) {
        std::cerr << "Failed to find geometry at " << geom_path << std::endl;
        return 1;
    }

    fvw::BladeGeometry geom = fvw::computeBladeGeometry(config.turbine, bladeDef);
    
    // 3. Define Output Paths
    std::string case_name = "NTNU_Baseline"; 
    // Usually this comes from config.perturbations or CLI. 
    // We assume the user wants to process the specific failed case.
    if (!config.perturbations.empty()) case_name = config.perturbations[0].name;

    std::string h5_file = "results/" + case_name + "/wake.h5";
    if (argc > 2) h5_file = argv[2];

    std::string output_csv = "results/" + case_name + "/probe_output_partial.csv";

    std::cout << "Processing H5 file: " << h5_file << std::endl;
    std::cout << "Writing CSV to: " << output_csv << std::endl;

    // 3b. Detect actual wake position (Hack for mismatched Hub Height)
    // Read t=0 wake and find average Z
    fvw::Wake tempWake(config.turbine.nBlades, config.turbine.nSegments, config.turbine.nSegments+1);
    try {
        fvw::read_wake_snapshot(tempWake, h5_file, 0, config.turbine);
        double sumZ = 0;
        int count = 0;
        for(int b=0; b<config.turbine.nBlades; ++b) {
            for(const auto& n : tempWake.getBladeWake(0, b).nodes) {
                sumZ += n.position.z;
                count++;
            }
        }
        if (count > 0) {
            double avgZ = sumZ / count;
            std::cout << "Detected Wake Center Z = " << avgZ << " m. (Ref Config HubHeight = " << config.turbine.hubHeight << " m)" << std::endl;
            if (std::abs(avgZ - config.turbine.hubHeight) > 1.0) {
                std::cout << "WARNING: Significant mismatch in Hub Height! Overriding config with detected height." << std::endl;
                config.turbine.hubHeight = avgZ;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Could not auto-detect wake position: " << e.what() << std::endl;
    }

    // 4. Run Probe Calculation
    fvw::runProbeCalculation(h5_file, output_csv, config.turbine, config.sim, geom);

    return 0;
}
