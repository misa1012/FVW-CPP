#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <cmath>
#include <filesystem>

#include "wake.h"
#include "postprocess.h"
#include "geometry.h"
#include "config.h"

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

    const int target_timestep = config.sim.timesteps - 1;

    // Load geometry (logic similar to others, reusable)
    std::string turbine_model = "NREL_5MW"; // should be from config
    std::string data_root = "data/" + turbine_model;
    std::string geometry_path = data_root + "/blade_geometry.csv";
    if (!fs::exists(geometry_path)) {
        data_root = "../data/" + turbine_model;
        geometry_path = data_root + "/blade_geometry.csv";
    }

    // Load geometry from file
    auto bladeDef = fvw::loadBladeDefinition(geometry_path);
    if (bladeDef.r.empty()) {
         std::cerr << "Fatal: Could not load blade geometry from " << geometry_path << "\n";
         return 1;
    }
    auto geom = fvw::computeBladeGeometry(config.turbine, bladeDef);
    std::cout << "Blade geometry computed.\n";

    // ====== 统一 1 m 网格（3D）======
    const double D = 2.0 * config.turbine.rTip;
    const double hub_height = config.turbine.hubHeight;
    const double h = 1.0; // 1 m

    const double x_start = -1.0 * D;
    const double x_end   =  3.0 * D;

    const double half_span = 1.5 * D / 2.0; // 0.75D
    const double y_min = -half_span;
    const double y_max =  half_span;

    const double z_min = hub_height - half_span;
    const double z_max = hub_height + half_span;

    std::cout << "Building 1 m uniform 3D grid...\n";
    Grid3D grid = build_uniform_grid_3d(x_start, x_end, y_min, y_max, z_min, z_max, h);
    std::cout << "Grid: Nx=" << grid.Nx << " Ny=" << grid.Ny << " Nz=" << grid.Nz
              << "  total=" << grid.points.size() << "\n";

    // ====== 读取 wake 并计算诱导速度 ======
    fvw::Wake wake_snapshot(config.turbine.nBlades, config.turbine.nSegments, config.turbine.nSegments + 1);
    fvw::read_wake_snapshot(wake_snapshot, h5_path.string(), target_timestep, config.turbine);

    if (wake_snapshot.bladeWakes.empty()) {
        std::cerr << "[ERR] Failed to load wake data from " << h5_path << "\n";
        return 1;
    }

    std::vector<fvw::Vec3> velocities(grid.points.size(), {0.0, 0.0, 0.0});
    auto t0 = std::chrono::high_resolution_clock::now();
    fvw::computeInducedVelocity(velocities, grid.points, wake_snapshot, target_timestep, config.turbine, geom, config.sim);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
    std::cout << "Induced velocity computed in " << sec << " s\n";

    // ====== 写 VTK ======
    write_field_to_vtk(out_vtk, grid.points, velocities, grid.Nx, grid.Ny, grid.Nz);
    std::cout << "Wrote VTK: " << out_vtk << "\n";

    return 0;
}
