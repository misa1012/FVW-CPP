#include "simulation/simulation_runner.h"
#include "io/cli_utils.h"
#include "io/logger.h"

#include <fstream>
#include <sstream>

namespace fvw {

namespace {
const char* to_string(VortexModelType t) {
    switch (t) {
        case VortexModelType::Constant: return "Constant";
        case VortexModelType::GammaDecay: return "GammaDecay";
    }
    return "Unknown";
}

const char* to_string(VortexCoreType t) {
    switch (t) {
        case VortexCoreType::VanGarrel: return "VanGarrel";
        case VortexCoreType::VanGarrelUnitConsistent: return "VanGarrelUnitConsistent";
        case VortexCoreType::ChordBasedCore: return "ChordBasedCore";
    }
    return "Unknown";
}

} // namespace

SimulationRunner::SimulationRunner(const PerturbationConfig& pc,
                                   const GlobalConfig& config,
                                   const std::string& rootOutput,
                                   const RunMetadata& runMeta)
    : m_pc(pc), m_globalConfig(config), m_rootOutput(rootOutput), m_runMeta(runMeta)
{
    m_turbineParams = m_globalConfig.turbine;
    m_simParams = m_globalConfig.sim;
}

void SimulationRunner::initialize() {
    // 1. Create case output directory
    m_caseOutDir = (std::filesystem::path(m_rootOutput) / m_pc.name).string();
    std::filesystem::create_directories(m_caseOutDir);
    m_h5Filepath = (std::filesystem::path(m_caseOutDir) / "wake.h5").string();

    // 2. Initialize Logger
    Logger::init(m_caseOutDir + "/simulation.log");
    Logger::info("Case", m_pc.name);
    Logger::info("Output Dir", m_caseOutDir);
    
    // 3. Resolve resource paths
    resolve_paths();

    // Print run summary to stdout (captured in .out by job script)
    print_run_summary_stdout();

    // 4. Apply perturbation
    apply_perturbation();

    // 5. Load resources (Geometry, Airfoils)
    load_resources();
    char buffer[100];
    snprintf(buffer, sizeof(buffer), "%.2f rad/s (%.1f rpm)", m_turbineParams.omega, m_turbineParams.omega * 60.0 / (2.0 * M_PI));
    Logger::info("Rotor Speed", buffer);
    Logger::info("TSR", std::to_string(m_turbineParams.tsr));
    
    int num_steps = m_simParams.timesteps - 1;
    snprintf(buffer, sizeof(buffer), "%.4f s (dt=%.6f s, %d steps)", m_simParams.totalTime, m_simParams.dt, num_steps);
    Logger::info("Simulation Time", buffer);
    
    snprintf(buffer, sizeof(buffer), "%d blades x %d segments", m_turbineParams.nBlades, m_turbineParams.nSegments);
    Logger::info("Discretization", buffer);

    // 6. Compute derived geometry
    m_geom = computeBladeGeometry(m_turbineParams, m_bladeDef);

    // 7. Allocate state containers
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
    computeBEM(*m_perf, m_geom, m_turbineParams, m_airfoils, m_simParams);

    InitializeWakeStructure(*m_wake, m_geom, *m_perf, m_turbineParams, *m_pos, m_simParams);
    
    // Write initial state
    writeWakeToHDF5(*m_wake, *m_perf, m_turbineParams, m_h5Filepath, 0);
    writeConfigToHDF5(m_geom, m_turbineParams, m_simParams, m_h5Filepath);
}

void SimulationRunner::run() {
    auto total_start = std::chrono::high_resolution_clock::now();
    
    cli::ProgressBar progress(m_simParams.timesteps);
    
    // Initial update
    progress.update(0, 0.0);

    for (int t = 1; t < m_simParams.timesteps; ++t)
    {
        double currentTime = t * m_simParams.dt;
        if (m_simParams.logStepTiming || m_simParams.logVerbose) {
            Logger::section_header(currentTime, t);
        }
        auto step_start = std::chrono::high_resolution_clock::now();

        AdvanceWakeStructure(*m_wake, m_turbineParams, *m_pos, m_simParams.dt, t);
        kuttaJoukowskiIteration(*m_wake, *m_perf, m_geom, *m_axes, *m_pos, *m_velBCS, m_airfoils, m_turbineParams, m_simParams);

        if (m_simParams.vortexModel == VortexModelType::GammaDecay)
        {
            ApplyGammaDecayAndRemoval(*m_wake, t, m_turbineParams, m_simParams);
        }

        UpdateWakeVelocities(*m_wake, m_turbineParams, t, m_geom, m_simParams);

        if (t % m_simParams.outputFrequency == 0 || t == m_simParams.timesteps - 1)
        {
            writeWakeToHDF5(*m_wake, *m_perf, m_turbineParams, m_h5Filepath, t);
        }

        auto step_end = std::chrono::high_resolution_clock::now();
        if (m_simParams.logStepTiming) {
            double step_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(step_end - step_start).count();
            std::ostringstream oss;
            oss << "Step " << t << " elapsed: " << std::fixed << std::setprecision(6) << step_seconds << " s";
            Logger::log("TIME", oss.str());
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
    m_totalSimSeconds = total_duration.count() / 1000.0;
    std::cout << cli::GREEN << "Done in " << m_totalSimSeconds << " s" << cli::RESET << std::endl;
}

void SimulationRunner::finalize() {
    Logger::info("Case '" + m_pc.name + "' completed.");
    if (m_simParams.logStepTiming) {
        std::ostringstream oss;
        oss << "Total simulation time: " << std::fixed << std::setprecision(6) << m_totalSimSeconds << " s";
        Logger::log("TIME", oss.str());
    }
    Logger::close();
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
    m_dataRoot = data_root;

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

void SimulationRunner::print_run_summary_stdout() const {
    std::cout << "=== Run Summary ===" << std::endl;
    if (!m_runMeta.exe_path.empty()) {
        std::cout << "Executable   : " << m_runMeta.exe_path << std::endl;
    }
    if (!m_runMeta.cli_args.empty()) {
        std::cout << "CLI Args     : " << m_runMeta.cli_args << std::endl;
    }
    if (!m_runMeta.config_path.empty()) {
        std::cout << "Config       : " << m_runMeta.config_path << std::endl;
    }
    std::cout << "Case         : " << m_pc.name << std::endl;
    std::cout << "Output Dir   : " << m_caseOutDir << std::endl;
    std::cout << "Data Root    : " << m_dataRoot << std::endl;
    std::cout << "Geometry     : " << m_geometryPath << std::endl;
    std::cout << "Model        : " << m_turbineParams.model << std::endl;
    std::cout << "Wind Speed   : " << m_turbineParams.windSpeed << " m/s" << std::endl;
    std::cout << "TSR          : " << m_turbineParams.tsr << std::endl;
    std::cout << "Omega        : " << m_turbineParams.omega << " rad/s" << std::endl;
    std::cout << "dt           : " << m_simParams.dt << " s" << std::endl;
    std::cout << "Timesteps    : " << m_simParams.timesteps << std::endl;
    std::cout << "Output Freq  : " << m_simParams.outputFrequency << std::endl;
    std::cout << "Core Type    : " << to_string(m_simParams.coreType) << std::endl;
    std::cout << "Vortex Model : " << to_string(m_simParams.vortexModel) << std::endl;
    std::cout << "===================" << std::endl;
}

} // namespace fvw
