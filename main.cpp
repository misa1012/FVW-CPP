#include "airfoil.h"
#include "geometry.h"
#include "position.h"
#include "velocity.h"
#include "performance.h"
#include "bem.h"
#include "wake.h"
#include <iostream>
#include <filesystem>

#include "validate.h"

int main(int argc, char *argv[])
{
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

    // 解析验证参数
    auto valParams = fvw::parseValidationArgs(argc, argv);

    // Section 1. Define parameters
    // Call geometry.h
    // Simulation parameters
    fvw::SimParams simParams;
    simParams.dt = 0.06;
    simParams.totalTime = 1.2;
    simParams.timesteps = static_cast<int>(simParams.totalTime / simParams.dt) + 1;
    // Turbine parameters
    fvw::TurbineParams turbineParams;
    turbineParams.windSpeed = 11.4;
    turbineParams.rho = 1.23;
    turbineParams.rHub = 1.5;
    turbineParams.rTip = 63.0;
    turbineParams.nBlades = 3;
    turbineParams.nSegments = 18;
    turbineParams.tsr = 7.0;
    turbineParams.omega = turbineParams.tsr * turbineParams.windSpeed / turbineParams.rTip;

    // Compute blade geometry
    auto geom = fvw::computeBladeGeometry(turbineParams);
    std::cout << "Blade geometry computed." << std::endl;

    // Section 2. Read in airfoil files
    // Call airfoil.h
    // Read airfoils
    auto airfoils = fvw::readAirfoils("../data/NREL_files/");
    std::cout << "Airfoil profiles read-in completed." << std::endl;

    // Section 3. Calculate the initial position
    // Call position.h
    // Compute positions
    fvw::PositionData pos(turbineParams.nBlades, simParams.timesteps,
                          turbineParams.nSegments + 1, turbineParams.nSegments);
    fvw::computePositions(pos, simParams, turbineParams, geom);
    std::cout << "Position computation completed." << std::endl;

    // 验证 geometry（如果需要）
    if (!valParams.section.empty() &&
        (valParams.section == "all" || valParams.section == "geometry"))
    {
        fvw::ValidationParams geomParams = valParams;
        geomParams.section = "geometry"; // 仅验证 geometry
        fvw::validate(geomParams, simParams, turbineParams, geom, pos);
    }

    // 验证 position（如果需要）
    if (!valParams.section.empty() &&
        (valParams.section == "all" || valParams.section == "position"))
    {
        fvw::ValidationParams posParams = valParams;
        posParams.section = "position"; // 仅验证 position
        fvw::validate(posParams, simParams, turbineParams, geom, pos);
    }

    // Call velocity.h
    // Compute velocities
    // Section 4. Compute velocities
    fvw::VelICS velICS(turbineParams.nBlades, simParams.timesteps, turbineParams.nSegments);
    fvw::computeVelICS(velICS, pos, simParams, turbineParams);

    fvw::VelBCS velBCS(turbineParams.nBlades, simParams.timesteps, turbineParams.nSegments);
    fvw::NodeAxes axes(turbineParams.nBlades, simParams.timesteps,
                       turbineParams.nSegments + 1, turbineParams.nSegments);
    fvw::computeVelBCS(velBCS, velICS, axes, pos, simParams, turbineParams);
    std::cout << "Velocity computation completed." << std::endl;

    // Initialize performance data
    fvw::PerformanceData perf(turbineParams.nBlades, simParams.timesteps, turbineParams.nSegments);
    std::cout << "Performance data is initialized." << std::endl;

    // Compute BEM
    fvw::computeBEM(perf, geom, turbineParams, airfoils);
    std::cout << "BEM computation completed." << std::endl;

    // Initialize the wake for the first timestep
    fvw::Wake wake(turbineParams.nBlades, turbineParams.nSegments, turbineParams.nSegments + 1);
    // t = 0
    initializeWake(wake, geom, perf, turbineParams, pos, simParams.dt);

    // 现在line的位置初始化好了，应该更新gamma，这是需要循环更新
    std::vector<fvw::VortexLine> linesLifting(wake.nBlades * wake.nShed);
    for (size_t i = 0; i < wake.lines[1].size(); ++i)
    {
        if (wake.lines[1][i].type == fvw::VortexLineType::Bound)
        {
            linesLifting.push_back(wake.lines[1][i]);
        }
    }
    kuttaJoukowskiIteration(linesLifting, wake.nodes[1], perf, geom, axes, turbineParams, pos, velBCS, airfoils);



    // // Update wake and compute induced velocity (示例循环)
    // for (int t = 1; t < simParams.timesteps; ++t)
    // {
    //     fvw::updateWake(wake, vel, simParams);
    //     fvw::computeInducedVelocity(vel, wake, geom, turbineParams, simParams);
    // fvw::computeVelICS(velICS, pos, simParams, turbineParams);
    // fvw::computeVelBCS(velBCS, velICS, axes, pos, simParams, turbineParams);
    //     // 重新计算迎角和 BEM
    //     fvw::computeAoAG(aoag, vel, turbineParams.nBlades, simParams.timesteps, turbineParams.nSegments);
    //     fvw::computeBEM(perf, aoag, geom, turbineParams, airfoils);
    // }
    // std::cout << "Wake computation completed." << std::endl;

    return 0;
}