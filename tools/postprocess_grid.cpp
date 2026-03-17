#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <cmath>
#include <filesystem>

#include "core/wake.h"
#include "io/postprocess.h"
#include "core/geometry.h"
#include "io/config.h"

namespace fs = std::filesystem;

// ---- 写 VTK ----
static void write_field_to_vtk(const std::string &filename,
                               const std::vector<fvw::Vec3> &grid_points,
                               const std::vector<fvw::Vec3> &velocities,
                               int Nx, int Ny, int Nz)
{
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open())
    {
        std::cerr << "[ERR] Unable to open VTK file " << filename << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Wake Velocity Field\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    vtkFile << "POINTS " << grid_points.size() << " float\n";
    for (const auto &p : grid_points)
        vtkFile << p.x << " " << p.y << " " << p.z << "\n";

    vtkFile << "POINT_DATA " << velocities.size() << "\n";
    vtkFile << "VECTORS Velocity float\n";
    for (const auto &v : velocities)
        vtkFile << v.x << " " << v.y << " " << v.z << "\n";

    vtkFile.close();
}

// ---- 生成 1m 统一 3D 网格 ----
struct Grid3D {
    std::vector<double> x, y, z;
    std::vector<fvw::Vec3> points;
    int Nx=0, Ny=0, Nz=0;
};

static Grid3D build_uniform_grid_3d(double x_start, double x_end,
                                    double y_min, double y_max,
                                    double z_min, double z_max,
                                    double h)
{
    Grid3D g;
    const double eps = 1e-9;

    for (double x = x_start; x <= x_end + eps; x += h) g.x.push_back(x);
    for (double y = y_min;   y <= y_max + eps; y += h) g.y.push_back(y);
    for (double z = z_min;   z <= z_max + eps; z += h) g.z.push_back(z);

    g.Nx = (int)g.x.size(); g.Ny = (int)g.y.size(); g.Nz = (int)g.z.size();

    g.points.reserve((size_t)g.Nx * g.Ny * g.Nz);
    for (int k = 0; k < g.Nz; ++k)
        for (int j = 0; j < g.Ny; ++j)
            for (int i = 0; i < g.Nx; ++i)
                g.points.push_back({g.x[i], g.y[j], g.z[k]});

    return g;
}

static Grid3D build_rotated_vertical_plane(double x_start, double x_end,
                                           double s_min, double s_max,
                                           double hub_height,
                                           double alpha_deg,
                                           double h)
{
    Grid3D g;
    const double eps = 1e-9;
    const double alpha = alpha_deg * M_PI / 180.0;
    const double sa = std::sin(alpha);
    const double ca = std::cos(alpha);

    for (double x = x_start; x <= x_end + eps; x += h) g.x.push_back(x);
    for (double s = s_min; s <= s_max + eps; s += h) g.y.push_back(s);
    g.z.push_back(0.0);

    g.Nx = (int)g.x.size();
    g.Ny = (int)g.y.size();
    g.Nz = 1;

    g.points.reserve((size_t)g.Nx * g.Ny);
    for (int j = 0; j < g.Ny; ++j) {
        const double s = g.y[j];
        const double y = s * sa;
        const double z = hub_height + s * ca;
        for (int i = 0; i < g.Nx; ++i) {
            g.points.push_back({g.x[i], y, z});
        }
    }

    return g;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " /path/to/wake.h5 [optional:/path/to/output.vtk] [optional:/path/to/config.json]\n";
        return 1;
    }

    const fs::path h5_path(argv[1]);
    if (!fs::exists(h5_path)) {
        std::cerr << "[ERR] wake.h5 not found: " << h5_path << "\n";
        return 1;
    }

    // 输出文件
    std::string out_vtk;
    if (argc >= 3) {
        out_vtk = argv[2];
    } else {
        out_vtk = (h5_path.parent_path() / "vel_grid_3d_final.vtk").string();
    }

    // Load Config
    std::string config_path = "config.json";
    if (argc >= 4) {
        config_path = argv[3];
    } else {
        // Try neighbor for config if not found
        if (!fs::exists(config_path) && fs::exists("../config.json")) config_path = "../config.json";
    }

    fvw::GlobalConfig config;
    try {
        config = fvw::ConfigLoader::load(config_path);
    } catch(const std::exception& e) {
        std::cerr << "Error loading config: " << e.what() << "\nUsing defaults.\n";
        // Fallback or exit? Should exit.
        return 1;
    }

    int target_timestep = config.sim.timesteps - 1;
    if (argc >= 5) {
        target_timestep = std::stoi(argv[4]);
    }
    
    double grid_res = 2.0; // Default to 2.0m to save space/time
    if (argc >= 6) {
        grid_res = std::stod(argv[5]);
    }

    int slice_mode = 0;
    if (argc >= 7) {
        slice_mode = std::stoi(argv[6]);
    }

    // Load geometry
    std::string turbine_model = config.turbine.model; // Use model from config
    if (turbine_model.empty()) turbine_model = "NREL_5MW"; // Fallback

    std::string data_root = "data/" + turbine_model;
    if (!fs::exists(data_root)) {
        if (fs::exists("../data/" + turbine_model)) {
            data_root = "../data/" + turbine_model;
        } else {
             std::cerr << "[WARN] Data directory for model " << turbine_model << " not found. Trying default locations.\n";
        }
    }
    
    std::string geometry_path = data_root + "/blade_geometry.csv";
    std::cout << "Loading geometry from: " << geometry_path << "\n";

    // Load geometry from file
    auto bladeDef = fvw::loadBladeDefinition(geometry_path);
    if (bladeDef.r.empty()) {
         std::cerr << "Fatal: Could not load blade geometry from " << geometry_path << "\n";
         return 1;
    }
    auto geom = fvw::computeBladeGeometry(config.turbine, bladeDef);
    std::cout << "Blade geometry computed.\n";

    // ====== 统一 1 m 网格（3D）======
    double x_start_factor = -1.0;
    if (argc >= 8) {
        x_start_factor = std::stod(argv[7]);
    }

    double x_end_factor = 6.0; // Default expanded to 6D
    if (argc >= 9) {
        x_end_factor = std::stod(argv[8]);
    }

    const double D = 2.0 * config.turbine.rTip;
    const double hub_height = config.turbine.hubHeight;
    const double h = grid_res;

    const double x_start = x_start_factor * D;
    const double x_end   = x_end_factor * D;
    
    std::cout << "[DEBUG] Grid Bounds Factors: " << x_start_factor << " to " << x_end_factor << " (D=" << D << ")\n";
    std::cout << "[DEBUG] Grid Physical X: " << x_start << " to " << x_end << " m\n";

    const double half_span = 3.0 * D / 2.0; // 1.5D (Covers -D to D with margin)
    
    // Default Volume
    double y_min = -half_span;
    double y_max =  half_span;
    double z_min = hub_height - half_span;
    double z_max = hub_height + half_span;
    
    double plane_angle_deg = 0.0;
    if (argc >= 10) {
        plane_angle_deg = std::stod(argv[9]);
    }

    if (slice_mode == 1) {
        // Horizontal Slice (Z fixed)
        z_min = hub_height;
        z_max = hub_height;
        std::cout << "Slice Mode 1: Horizontal Plane (Z = " << hub_height << ")\n";
    } else if (slice_mode == 2) {
        // Vertical Slice (Y fixed)
        y_min = 0.0;
        y_max = 0.0;
        std::cout << "Slice Mode 2: Vertical Plane (Y = 0.0)\n";
    } else if (slice_mode == 3) {
        std::cout << "Slice Mode 3: Rotated Vertical Plane, alpha = " << plane_angle_deg << " deg\n";
    } else {
        std::cout << "Slice Mode 0: Full 3D Volume\n";
    }

    std::cout << "Building grid...\n";
    Grid3D grid;
    if (slice_mode == 3) {
        grid = build_rotated_vertical_plane(x_start, x_end, -half_span, half_span, hub_height, plane_angle_deg, h);
    } else {
        grid = build_uniform_grid_3d(x_start, x_end, y_min, y_max, z_min, z_max, h);
    }
    std::cout << "Grid: Nx=" << grid.Nx << " Ny=" << grid.Ny << " Nz=" << grid.Nz
              << "  total=" << grid.points.size() << "\n";

    // ====== 读取 wake 并计算诱导速度 ======
    fvw::Wake wake_snapshot(config.turbine.nBlades, config.turbine.nSegments, config.turbine.nSegments + 1);
    fvw::read_wake_snapshot(wake_snapshot, h5_path.string(), target_timestep, config.turbine);

    if (wake_snapshot.bladeWakes.empty()) {
        std::cerr << "[ERR] Wake snapshot has no blade wakes (unexpected structure).\n";
        return 1;
    }

    size_t totalNodes = 0;
    size_t totalLines = 0;
    double totalGamma = 0.0;
    for (int b = 0; b < config.turbine.nBlades; ++b) {
        const auto& bw = wake_snapshot.getBladeWake(target_timestep, b);
        totalNodes += bw.nodes.size();
        totalLines += bw.lines.size();
        for (const auto& line : bw.lines) {
            totalGamma += std::abs(line.gamma);
        }
    }
    std::cout << "[DEBUG] Loaded Wake: " << totalNodes << " nodes, " << totalLines << " lines. Total |Gamma|: " << totalGamma << "\n";

    if (totalNodes == 0 || totalLines == 0) {
        std::cerr << "[WARN] Wake appears empty! Induced velocity will be ZERO.\n";
    }

    std::vector<fvw::Vec3> velocities(grid.points.size(), {0.0, 0.0, 0.0});
    auto t0 = std::chrono::high_resolution_clock::now();
    fvw::computeInducedVelocity(velocities, grid.points, wake_snapshot, target_timestep, geom, config.sim);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
    std::cout << "Induced velocity computed in " << sec << " s\n";

    // Add free stream velocity
    double U_inf = config.turbine.windSpeed;
    std::cout << "Adding free stream velocity: " << U_inf << " m/s\n";
    for (auto &v : velocities) {
        v.x += U_inf;
    }

    // ====== 写 VTK ======
    write_field_to_vtk(out_vtk, grid.points, velocities, grid.Nx, grid.Ny, grid.Nz);
    std::cout << "Wrote VTK: " << out_vtk << "\n";

    return 0;
}
