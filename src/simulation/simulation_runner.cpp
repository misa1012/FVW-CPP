#include "simulation/simulation_runner.h"
#include "io/cli_utils.h"
#include "io/logger.h"

#include <array>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unistd.h>

namespace fvw {

namespace {
std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '\"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    std::ostringstream oss;
                    oss << "\\u" << std::hex << std::uppercase << std::setw(4)
                        << std::setfill('0') << static_cast<int>(static_cast<unsigned char>(c));
                    out += oss.str();
                } else {
                    out += c;
                }
        }
    }
    return out;
}

std::string now_iso8601() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
    localtime_r(&t, &tm);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S%z");
    return oss.str();
}

std::string get_hostname() {
    std::array<char, 256> buf{};
    if (gethostname(buf.data(), buf.size()) == 0) {
        return std::string(buf.data());
    }
    return "unknown";
}

std::string read_file_trim(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    std::string s = ss.str();
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

std::string read_git_commit(const std::filesystem::path& project_root) {
    std::filesystem::path head = project_root / ".git" / "HEAD";
    if (!std::filesystem::exists(head)) return {};
    std::string head_content = read_file_trim(head);
    const std::string ref_prefix = "ref: ";
    if (head_content.rfind(ref_prefix, 0) == 0) {
        std::filesystem::path ref_path = project_root / ".git" / head_content.substr(ref_prefix.size());
        return read_file_trim(ref_path);
    }
    return head_content;
}

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
        case VortexCoreType::ChordBasedCore: return "ChordBasedCore";
    }
    return "Unknown";
}

const char* to_string(SegmentDistribution t) {
    switch (t) {
        case SegmentDistribution::Linear: return "Linear";
        case SegmentDistribution::Cosine: return "Cosine";
    }
    return "Unknown";
}

const char* to_string(PerturbationType t) {
    switch (t) {
        case PerturbationType::None: return "None";
        case PerturbationType::CollectivePitch: return "CollectivePitch";
        case PerturbationType::AsymmetricStaticPitch: return "AsymmetricStaticPitch";
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

    // Save config and manifest early for traceability
    write_run_manifest();

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

void SimulationRunner::write_run_manifest() {
    // Build manifest JSON (includes raw config)
    std::filesystem::path root_path = std::filesystem::path(m_rootOutput);
    std::filesystem::path project_root = (root_path.filename() == "results") ? root_path.parent_path() : root_path;

    std::string git_commit = read_git_commit(project_root);
    std::string host = get_hostname();
    std::string timestamp = now_iso8601();
    std::string config_abs;
    if (!m_runMeta.config_path.empty()) {
        std::error_code ec;
        config_abs = std::filesystem::absolute(m_runMeta.config_path, ec).string();
    }

    std::filesystem::path manifest_path = std::filesystem::path(m_caseOutDir) / "run_manifest.json";
    std::ofstream out(manifest_path);
    if (!out) return;

    out << "{\n";
    out << "  \"timestamp\": \"" << json_escape(timestamp) << "\",\n";
    out << "  \"hostname\": \"" << json_escape(host) << "\",\n";
    out << "  \"executable\": \"" << json_escape(m_runMeta.exe_path) << "\",\n";
    out << "  \"cli_args\": \"" << json_escape(m_runMeta.cli_args) << "\",\n";
    out << "  \"config\": {\n";
    out << "    \"path\": \"" << json_escape(m_runMeta.config_path) << "\",\n";
    out << "    \"path_abs\": \"" << json_escape(config_abs) << "\",\n";
    out << "    \"raw_json\": \"" << json_escape(m_runMeta.config_text) << "\"\n";
    out << "  },\n";
    out << "  \"git\": {\n";
    out << "    \"commit\": \"" << json_escape(git_commit) << "\"\n";
    out << "  },\n";
    out << "  \"output\": {\n";
    out << "    \"root\": \"" << json_escape(m_rootOutput) << "\",\n";
    out << "    \"case_dir\": \"" << json_escape(m_caseOutDir) << "\",\n";
    out << "    \"wake_h5\": \"" << json_escape(m_h5Filepath) << "\",\n";
    out << "    \"log\": \"simulation.log\"\n";
    out << "  },\n";
    out << "  \"data\": {\n";
    out << "    \"data_root\": \"" << json_escape(m_dataRoot) << "\",\n";
    out << "    \"geometry_path\": \"" << json_escape(m_geometryPath) << "\"\n";
    out << "  },\n";
    out << "  \"case\": {\n";
    out << "    \"name\": \"" << json_escape(m_pc.name) << "\",\n";
    out << "    \"perturbation_type\": \"" << json_escape(to_string(m_pc.type)) << "\",\n";
    out << "    \"perturbation_amplitude_deg\": " << m_pc.amplitude_deg << ",\n";
    out << "    \"perturbation_freq_factor\": " << m_pc.freqFactor << "\n";
    out << "  },\n";
    out << "  \"turbine\": {\n";
    out << "    \"model\": \"" << json_escape(m_turbineParams.model) << "\",\n";
    out << "    \"wind_speed\": " << m_turbineParams.windSpeed << ",\n";
    out << "    \"rho\": " << m_turbineParams.rho << ",\n";
    out << "    \"r_hub\": " << m_turbineParams.rHub << ",\n";
    out << "    \"r_tip\": " << m_turbineParams.rTip << ",\n";
    out << "    \"hub_height\": " << m_turbineParams.hubHeight << ",\n";
    out << "    \"n_blades\": " << m_turbineParams.nBlades << ",\n";
    out << "    \"n_segments\": " << m_turbineParams.nSegments << ",\n";
    out << "    \"segment_distribution\": \"" << json_escape(to_string(m_turbineParams.segmentDistribution)) << "\",\n";
    out << "    \"tsr\": " << m_turbineParams.tsr << ",\n";
    out << "    \"omega\": " << m_turbineParams.omega << "\n";
    out << "  },\n";
    out << "  \"simulation\": {\n";
    out << "    \"dt\": " << m_simParams.dt << ",\n";
    out << "    \"total_time\": " << m_simParams.totalTime << ",\n";
    out << "    \"timesteps\": " << m_simParams.timesteps << ",\n";
    out << "    \"output_frequency\": " << m_simParams.outputFrequency << ",\n";
    out << "    \"cutoff_param\": " << m_simParams.cutoffParam << ",\n";
    out << "    \"vortex_model\": \"" << json_escape(to_string(m_simParams.vortexModel)) << "\",\n";
    out << "    \"core_type\": \"" << json_escape(to_string(m_simParams.coreType)) << "\"\n";
    out << "  }\n";
    out << "}\n";
}

} // namespace fvw
