#include "simulation/simulation_runner.h"
#include "io/cli_utils.h"

namespace fvw {

SimulationRunner::SimulationRunner(const PerturbationConfig& pc, const GlobalConfig& config, const std::string& rootOutput)
    : m_pc(pc), m_globalConfig(config), m_rootOutput(rootOutput)
{
    m_turbineParams = m_globalConfig.turbine;
    m_simParams = m_globalConfig.sim;
}

void SimulationRunner::initialize() {
    // 1. Resolve output paths
    m_caseOutDir = reset_case_output(m_rootOutput, m_pc.name, true);
    m_h5Filepath = (std::filesystem::path(m_caseOutDir) / "wake.h5").string();
    
    cli::print_header("Case: " + m_pc.name);
    cli::print_info("Output Dir", m_caseOutDir);
    
    // 2. Resolve resource paths
    resolve_paths();

    // 3. Apply perturbation
    apply_perturbation();

    // 4. Load resources (Geometry, Airfoils)
    load_resources();
    cli::print_info("Turbine", m_turbineParams.model);
    cli::print_info("Wind Speed", std::to_string(m_turbineParams.windSpeed) + " m/s");

    // 5. Compute derived geometry
    m_geom = computeBladeGeometry(m_turbineParams, m_bladeDef);

    // 6. Allocate state containers
    m_pos = std::make_unique<PositionData>(m_turbineParams.nBlades, m_simParams.timesteps,
                          m_turbineParams.nSegments + 1, m_turbineParams.nSegments);
    
    m_velICS = std::make_unique<VelICS>(m_turbineParams.nBlades, m_simParams.timesteps, m_turbineParams.nSegments);
    m_velBCS = std::make_unique<VelBCS>(m_turbineParams.nBlades, m_simParams.timesteps, m_turbineParams.nSegments);
    m_axes = std::make_unique<NodeAxes>(m_turbineParams.nBlades, m_simParams.timesteps,
                       m_turbineParams.nSegments + 1, m_turbineParams.nSegments);
    m_perf = std::make_unique<PerformanceData>(m_turbineParams.nBlades, m_simParams.timesteps, m_turbineParams.nSegments);
    m_wake = std::make_unique<Wake>(m_turbineParams.nBlades, m_turbineParams.nSegments, m_turbineParams.nSegments + 1);

    // 7. Initial Computation (t=0)
    computePositions(*m_pos, m_simParams, m_turbineParams, m_geom);
    computeVelICS(*m_velICS, *m_pos, m_simParams, m_turbineParams);
    computeVelBCS(*m_velBCS, *m_velICS, *m_axes, *m_pos, m_simParams, m_turbineParams);
    computeBEM(*m_perf, m_geom, m_turbineParams, m_airfoils);

    InitializeWakeStructure(*m_wake, m_geom, *m_perf, m_turbineParams, *m_pos, m_simParams);
    
    // Write initial state
    writeWakeToHDF5(*m_wake, *m_pos, *m_perf, *m_velICS, *m_velBCS, m_turbineParams, m_h5Filepath, 0);
    writeConfigToHDF5(m_geom, m_turbineParams, m_simParams, m_h5Filepath);
}

void SimulationRunner::run() {
    auto total_start = std::chrono::high_resolution_clock::now();
    
    cli::ProgressBar progress(m_simParams.timesteps);
    
    // Initial update
    progress.update(0, 0.0);

    for (int t = 1; t < m_simParams.timesteps; ++t)
    {
        AdvanceWakeStructure(*m_wake, m_geom, m_turbineParams, *m_pos, m_simParams.dt, t);
        kuttaJoukowskiIteration(*m_wake, *m_perf, m_geom, *m_axes, m_turbineParams, *m_pos, *m_velBCS, m_airfoils, m_simParams);

        if (m_simParams.vortexModel == VortexModelType::GammaDecay)
        {
            ApplyGammaDecayAndRemoval(*m_wake, t, m_turbineParams, m_simParams);
        }

        UpdateWakeVelocities(*m_wake, m_turbineParams, t, m_geom, m_simParams);

        if (t % m_simParams.outputFrequency == 0 || t == m_simParams.timesteps - 1)
        {
            writeWakeToHDF5(*m_wake, *m_pos, *m_perf, *m_velICS, *m_velBCS, m_turbineParams, m_h5Filepath, t);
        }

        // Update progress bar every 10 steps or if done
        if (t % 10 == 0 || t == m_simParams.timesteps - 1) {
            auto current_now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_now - total_start).count() / 1000.0;
            progress.update(t, elapsed);
        }
    }
    progress.finish();

    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - total_start);
    std::cout << cli::GREEN << "Done in " << total_duration.count() / 1000.0 << " s" << cli::RESET << std::endl;
}

void SimulationRunner::finalize(bool projectToGrid, bool computeProbes) {
    if (projectToGrid)
    {
        std::cout << "\n--- Projecting final wake to Eulerian grid ---" << std::endl;
        SimParams sim_copy = m_simParams;
        sim_copy.cutoffParam = 1.0; 
        projectWakeToEulerianGrid(*m_wake, m_turbineParams, sim_copy, m_geom, m_simParams.timesteps - 1, m_caseOutDir);
    }

    if (computeProbes)
    {
        std::string probe_csv_filepath = (std::filesystem::path(m_caseOutDir) / "probe_output.csv").string();
        runProbeCalculation(m_h5Filepath, probe_csv_filepath, m_turbineParams, m_simParams, m_geom);
    }
    
    std::cout << "\n[END] Case '" << m_pc.name << "' completed." << std::endl;
}

// --- Helpers ---

void SimulationRunner::resolve_paths() {
    // Use model name from configuration
    std::string turbine_model = m_turbineParams.model;
    if (turbine_model.empty()) turbine_model = "NREL_5MW"; // Fallback 
    std::string data_root = "data/" + turbine_model;
    
    // Check local data first, then parent (for build/ execution)
    if (!std::filesystem::exists(data_root)) {
        if (std::filesystem::exists("../data/" + turbine_model)) {
            data_root = "../data/" + turbine_model;
        }
    }

    m_geometryPath = data_root + "/blade_geometry.csv";
    
    if (!std::filesystem::exists(m_geometryPath)) {
         std::cerr << "Fatal: Geometry file not found: " << m_geometryPath << "\n";
         exit(1);
    }
}

void SimulationRunner::load_resources() {
    // Load Geometry definition
    static BladeDefinition cachedBladeDef; // Optimization: load once if same file? 
    // Careful with static if multiple runners run in parallel or different models. 
    // For now, let's load every time to be safe, or check if empty.
    
    if (cachedBladeDef.r.empty()) {
         cachedBladeDef = loadBladeDefinition(m_geometryPath);
         if (cachedBladeDef.r.empty()) {
             std::cerr << "Fatal: Could not load blade geometry from " << m_geometryPath << "\n";
             exit(1);
         }
    }
    m_bladeDef = cachedBladeDef;

    // Load Airfoils
    // Path logic for airfoils needs to be robust too
    std::filesystem::path geomP(m_geometryPath);
    std::string data_root = geomP.parent_path().string();
    std::string airfoil_list_path = data_root + "/airfoil_list.txt";
    std::string airfoil_dir = data_root + "/airfoils/";
    
    m_airfoils = readAirfoils(airfoil_dir, airfoil_list_path);
    if (m_airfoils.empty()) {
        std::cerr << "Fatal: Could not load airfoils from " << airfoil_list_path << "\n";
        exit(1);
    }
}

void SimulationRunner::apply_perturbation() {
    m_simParams.perturbation.type = m_pc.type;
    m_simParams.perturbation.amplitude_deg = m_pc.amplitude_deg;

    if (m_pc.freqFactor > 0.0)
    {
        const double rotational_freq_hz = m_turbineParams.omega / (2.0 * M_PI);
        m_simParams.perturbation.frequency_hz = m_pc.freqFactor * rotational_freq_hz;
    }
    else
    {
        m_simParams.perturbation.frequency_hz = 0.0;
    }
}

std::string SimulationRunner::reset_case_output(const std::string &root, const std::string &case_name, bool clean)
{
    std::filesystem::path p(root);
    p /= case_name;
    if (case_name.empty() || case_name == "." || case_name == "..") return root; 
    
    if (clean && std::filesystem::exists(p))
    {
        std::error_code ec; 
        std::filesystem::remove_all(p, ec);
    }
    std::filesystem::create_directories(p);
    return p.string();
}

} // namespace fvw
