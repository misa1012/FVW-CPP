#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "core/geometry.h"
#include "core/wake.h"
#include "io/config.h"
#include "io/json_utils.h"
#include "io/postprocess.h"

namespace fs = std::filesystem;

namespace {
struct SamplePoint {
    std::string name;
    fvw::Vec3 pos;
};

std::string read_string(const fvw::json::Value& obj, const std::string& key, const std::string& def = "") {
    if (obj.contains(key)) return obj[key].as_string();
    return def;
}

int read_int(const fvw::json::Value& obj, const std::string& key, int def) {
    if (obj.contains(key)) return obj[key].as_int();
    return def;
}

double read_double(const fvw::json::Value& obj, const std::string& key, double def) {
    if (obj.contains(key)) return obj[key].as_double();
    return def;
}

fvw::Vec3 parse_vec3(const fvw::json::Value& arr) {
    const auto& a = arr.as_array();
    if (a.size() != 3) throw std::runtime_error("Expected 3-element array for vector");
    return {a[0].as_double(), a[1].as_double(), a[2].as_double()};
}

std::string resolve_geometry_path(const fvw::GlobalConfig& config) {
    std::string model = config.turbine.model;
    if (model.empty()) model = "NREL_5MW";
    std::string data_root = "data/" + model;
    std::string geom_path = data_root + "/blade_geometry.csv";
    if (!fs::exists(geom_path)) {
        data_root = "../data/" + model;
        geom_path = data_root + "/blade_geometry.csv";
    }
    return geom_path;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " sampling.json\n";
        return 1;
    }

    std::string sampling_path = argv[1];
    fvw::json::Value root = fvw::json::parse_file(sampling_path);

    std::string wake_h5 = read_string(root, "wake_h5");
    if (wake_h5.empty()) {
        std::cerr << "Missing required field: wake_h5\n";
        return 1;
    }

    std::string config_path = read_string(root, "config", "config.json");
    if (!fs::exists(config_path) && fs::exists("../" + config_path)) {
        config_path = "../" + config_path;
    }

    fvw::GlobalConfig config;
    try {
        config = fvw::ConfigLoader::load(config_path);
    } catch (const std::exception& e) {
        std::cerr << "Error loading config: " << e.what() << "\n";
        return 1;
    }

    std::string output_csv = read_string(root, "output_csv");
    if (output_csv.empty()) {
        fs::path h5p(wake_h5);
        output_csv = (h5p.parent_path() / "sampling_output.csv").string();
    }

    std::string units = read_string(root, "units", "m");
    double unit_scale = 1.0;
    if (units == "D" || units == "d") {
        unit_scale = 2.0 * config.turbine.rTip;
    }

    std::vector<SamplePoint> points;

    if (root.contains("points")) {
        const auto& arr = root["points"].as_array();
        int idx = 0;
        for (const auto& p : arr) {
            std::string name = read_string(p, "name");
            if (name.empty()) {
                name = "p" + std::to_string(idx);
            }
            fvw::Vec3 v{p["x"].as_double(), p["y"].as_double(), p["z"].as_double()};
            v.x *= unit_scale; v.y *= unit_scale; v.z *= unit_scale;
            points.push_back({name, v});
            idx++;
        }
    }

    if (root.contains("lines")) {
        const auto& arr = root["lines"].as_array();
        int line_idx = 0;
        for (const auto& l : arr) {
            std::string name = read_string(l, "name");
            if (name.empty()) {
                name = "line" + std::to_string(line_idx);
            }
            fvw::Vec3 start = parse_vec3(l["start"]);
            fvw::Vec3 end = parse_vec3(l["end"]);
            int n = read_int(l, "n", 2);
            if (n < 2) n = 2;
            start.x *= unit_scale; start.y *= unit_scale; start.z *= unit_scale;
            end.x *= unit_scale; end.y *= unit_scale; end.z *= unit_scale;

            for (int i = 0; i < n; ++i) {
                double t = static_cast<double>(i) / static_cast<double>(n - 1);
                fvw::Vec3 v{
                    start.x + t * (end.x - start.x),
                    start.y + t * (end.y - start.y),
                    start.z + t * (end.z - start.z)
                };
                points.push_back({name + "_" + std::to_string(i), v});
            }
            line_idx++;
        }
    }

    if (points.empty()) {
        std::cerr << "No sampling points defined (points/lines empty).\n";
        return 1;
    }

    // Geometry for induced velocity calculation
    std::string geom_path = resolve_geometry_path(config);
    auto bladeDef = fvw::loadBladeDefinition(geom_path);
    if (bladeDef.r.empty()) {
        std::cerr << "Fatal: Could not load blade geometry from " << geom_path << "\n";
        return 1;
    }
    auto geom = fvw::computeBladeGeometry(config.turbine, bladeDef);

    std::string mode = read_string(root, "mode", "instant");
    int total_steps = config.sim.timesteps;
    int step = total_steps - 1;

    std::vector<int> steps;
    if (mode == "average") {
        int start_step = total_steps - 50;
        if (start_step < 0) start_step = 0;
        int end_step = total_steps - 1;
        int stride = 1;

        if (root.contains("average")) {
            const auto& avg = root["average"];
            start_step = read_int(avg, "start_step", start_step);
            end_step = read_int(avg, "end_step", end_step);
            stride = read_int(avg, "stride", stride);
        }
        if (stride < 1) stride = 1;
        if (end_step < start_step) end_step = start_step;
        for (int t = start_step; t <= end_step; t += stride) {
            steps.push_back(t);
        }
    } else {
        if (root.contains("instant")) {
            step = read_int(root["instant"], "step", step);
        }
        steps.push_back(step);
    }

    fvw::Wake wake(config.turbine.nBlades, config.turbine.nSegments, config.turbine.nSegments + 1);
    std::vector<fvw::Vec3> sample_positions;
    sample_positions.reserve(points.size());
    for (const auto& p : points) sample_positions.push_back(p.pos);

    std::ofstream out(output_csv);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open output file: " << output_csv << "\n";
        return 1;
    }

    if (mode == "average") {
        std::vector<fvw::Vec3> sum(points.size(), {0.0, 0.0, 0.0});
        int n_samples = 0;
        for (int t : steps) {
            fvw::read_wake_snapshot(wake, wake_h5, t, config.turbine);
            std::vector<fvw::Vec3> induced;
            fvw::computeInducedVelocity(induced, sample_positions, wake, t, geom, config.sim);
            for (size_t i = 0; i < induced.size(); ++i) {
                sum[i].x += induced[i].x;
                sum[i].y += induced[i].y;
                sum[i].z += induced[i].z;
            }
            n_samples++;
        }

        out << "name,x,y,z,u_mean,v_mean,w_mean,samples\n";
        for (size_t i = 0; i < points.size(); ++i) {
            const auto& p = points[i];
            const auto& s = sum[i];
            out << p.name << ","
                << p.pos.x << "," << p.pos.y << "," << p.pos.z << ","
                << (n_samples ? s.x / n_samples : 0.0) << ","
                << (n_samples ? s.y / n_samples : 0.0) << ","
                << (n_samples ? s.z / n_samples : 0.0) << ","
                << n_samples << "\n";
        }
    } else {
        out << "step,time,name,x,y,z,u,v,w\n";
        for (int t : steps) {
            fvw::read_wake_snapshot(wake, wake_h5, t, config.turbine);
            std::vector<fvw::Vec3> induced;
            fvw::computeInducedVelocity(induced, sample_positions, wake, t, geom, config.sim);
            double time = t * config.sim.dt;
            for (size_t i = 0; i < points.size(); ++i) {
                const auto& p = points[i];
                const auto& v = induced[i];
                out << t << "," << time << ","
                    << p.name << ","
                    << p.pos.x << "," << p.pos.y << "," << p.pos.z << ","
                    << v.x << "," << v.y << "," << v.z << "\n";
            }
        }
    }

    std::cout << "Sampling output written to: " << output_csv << "\n";
    return 0;
}
