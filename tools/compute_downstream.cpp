#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>

#include "core/wake.h"
#include "io/postprocess.h"
#include "core/geometry.h"

namespace fs = std::filesystem;

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>

#include "core/wake.h"
#include "io/postprocess.h"
#include "core/geometry.h"
#include "io/config.h"

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    std::string config_path = "config.json";
    if (argc > 1) config_path = argv[1];

    if (!fs::exists(config_path)) {
       if (fs::exists("../" + config_path)) config_path = "../" + config_path;
       else {
           std::cerr << "[ERR] Config file not found: " << config_path << "\n";
           return 1;
       }
    }
    std::cout << "Loading config from " << config_path << "...\n";
    fvw::GlobalConfig config = fvw::ConfigLoader::load(config_path);

    // Use relative path for portability, relative to wherever we are running (usually project root or build)
    static std::string BASE_DIR = "./results";

    // Infer case folders from config perturbations
    std::vector<std::string> CASE_FOLDERS;
    if (config.perturbations.empty()) {
        CASE_FOLDERS.push_back("default_baseline"); // or whatever logic fvw_main uses
    } else {
        for(const auto& p : config.perturbations) {
            CASE_FOLDERS.push_back(p.name);
        }
    }

    // Downstream positions: 0.10D ~ 3.00D
    std::vector<double> XD_LIST;
    for (int i = 10; i <= 300; ++i) XD_LIST.push_back(i / 100.0);

    // 网格与物理常量 (From Config)
    const double R     = config.turbine.rTip;
    const double D     = 2.0 * R;
    const double HUB_H = config.turbine.hubHeight;
    const double U_INF = config.turbine.windSpeed;
    const double DT    = config.sim.dt; 
    
    // Time window for averaging (Last 50 steps from totalTime?)
    // Hardcoded logic in original file was: STEP_START=2500 (approx end), N_SNAPS=50.
    // Let's adapt to config.sim.timesteps.
    // Take last 50 steps.
    int total_steps = config.sim.timesteps;
    const int N_SNAPS = 50;
    const int STEP_STRIDE = 10; // maybe stride?
    // Ensure we don't go out of bounds
    int start_step = total_steps - (N_SNAPS * STEP_STRIDE);
    if (start_step < 0) start_step = total_steps - N_SNAPS; // fallback stride 1
    if (start_step < 0) start_step = 0; 

    // Override with original logic if totalTime is large enough?
    // The original code had static 2500 start. 
    // Let's stick to "Last portion of simulation" logic which is more generic.
    std::cout << "Analysing last " << N_SNAPS << " snapshots...\n";

    // 竖剖面 z 方向离散
    const double Z_MIN = -10.0;
    const double Z_MAX = HUB_H + R + 30.0; // Dynamic Range
    const double Z_DZ  = 1.0;

    // ROI 定义 ( Normalized by D )
    auto in_top_roi = [](double z_over_D) { return (z_over_D > 0.45 && z_over_D < 0.60); };
    auto in_bottom_roi = [](double z_over_D) { return (z_over_D > -0.60 && z_over_D < -0.45); };

    // Load geometry (Use config model name if we had it, or hardcode NREL/fallback logic)
    std::string turbine_model = "NREL_5MW"; // Should be in config
    if (config.turbine.rTip > 60.0) turbine_model = "NREL_5MW"; // heuristic
    // Use the same logic as fvw_main
    std::string data_root = "data/" + turbine_model;
    std::string geometry_path = data_root + "/blade_geometry.csv";
    if (!fs::exists(geometry_path)) {
        data_root = "../data/" + turbine_model;
        geometry_path = data_root + "/blade_geometry.csv";
    }

    auto bladeDef = fvw::loadBladeDefinition(geometry_path);
    if (bladeDef.r.empty()) {
        std::cerr << "Fatal: Could not load blade geometry from " << geometry_path << "\n";
        return 1;
    }
    
    auto geom = fvw::computeBladeGeometry(config.turbine, bladeDef);

    // Z-line
    std::vector<double> z_line;
    for (double z = Z_MIN; z <= Z_MAX + 1e-9; z += Z_DZ) z_line.push_back(z);

    for (const auto& folder : CASE_FOLDERS) {
        fs::path case_dir = fs::path(BASE_DIR) / folder;
        std::string h5 = (case_dir / "wake.h5").string();
        if (!fs::exists(h5)) {
            std::cerr << "[WARN] " << h5 << " not found, skipping.\n";
            continue;
        }

        fs::path out_csv = case_dir / "roi_points_y001.csv";
        std::ofstream ofs(out_csv);
        if (!ofs.is_open()) {
            std::cerr << "[ERR] Cannot write " << out_csv << "\n";
            continue;
        }

        ofs << "step,t,x_over_D,region,z_over_D,z_absolute,"
               "U_ind,V_ind,W_ind,U_total,u_perp\n";

        std::cout << "\n=== Case: " << folder << " ===\n";

        for (double xD : XD_LIST) {
            const double x_loc = xD * D;
            // std::cout << "  x/D = " << std::fixed << std::setprecision(2) << xD << " ...\r" << std::flush;

            std::vector<fvw::Vec3> grid;
            grid.reserve(z_line.size());
            for (double z : z_line) grid.push_back({x_loc, 0.0, z});

            for (int i = 0; i < N_SNAPS; ++i) {
                const int    step = start_step + i * STEP_STRIDE;
                if (step >= config.sim.timesteps) break;
                const double t    = step * DT;

                fvw::Wake wake(config.turbine.nBlades, config.turbine.nSegments, config.turbine.nSegments + 1);
                fvw::read_wake_snapshot(wake, h5, step, config.turbine);

                std::vector<fvw::Vec3> vel(grid.size());
                fvw::computeInducedVelocity(vel, grid, wake, step, config.turbine, geom, config.sim);

                for (size_t k = 0; k < grid.size(); ++k) {
                    const double z_abs     = grid[k].z;
                    const double z_over_D  = (z_abs - HUB_H) / D;
                    
                    bool is_top = in_top_roi(z_over_D);
                    bool is_bot = in_bottom_roi(z_over_D);

                    if (!is_top && !is_bot) continue;

                    const char* region = is_top ? "top" : "bottom";
                    const double Uind = vel[k].x;
                    const double V    = vel[k].y;
                    const double W    = vel[k].z;
                    const double Ut   = U_INF + Uind;
                    const double up   = std::sqrt(V*V + W*W);

                    ofs << step << ","
                        << std::fixed << std::setprecision(6) << t << ","
                        << std::fixed << std::setprecision(2) << xD << ","
                        << region << ","
                        << std::fixed << std::setprecision(6) << z_over_D << "," << z_abs << ","
                        << Uind << "," << V << "," << W << ","
                        << Ut << "," << up << "\n";
                }
            }
        }
        std::cout << "Wrote: " << out_csv << "\n";
        ofs.close();
    }

    std::cout << "\nAll Done.\n";
    return 0;
}
