#include "airfoil.h"
#include "geometry.h"
#include "position.h"
#include <iostream>
#include <filesystem>

int main() {
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

    // Simulation parameters
    fvw::SimParams simParams;
    simParams.dt = 0.06;
    simParams.totalTime = 12.0;
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

    // Read airfoils
    std::vector<std::string> airfoilNames = {"Cylinder1", "Cylinder2", "DU40_A17", "DU35_A17",
                                            "DU30_A17", "DU25_A17", "DU21_A17", "NACA64_A17"};
    auto airfoils = fvw::readAirfoils(airfoilNames);
    std::cout << "Airfoil profiles read-in completed." << std::endl;

    // Compute blade geometry
    auto geom = fvw::computeBladeGeometry(turbineParams);

    // Compute positions
    fvw::PositionData pos(turbineParams.nBlades, simParams.timesteps,
                          turbineParams.nSegments + 1, turbineParams.nSegments);
    fvw::computePositions(pos, simParams, turbineParams, geom);
    std::cout << "Position computation completed." << std::endl;

    // TODO: Add velocity, BEM, and wake modules

    return 0;
}