#ifndef FVW_AIRFOIL_H
#define FVW_AIRFOIL_H

#include <string>
#include <vector>

namespace fvw
{

    struct AirfoilData
    {
        double stallAoA;
        double cn0AoA;
        double lift0Cn;
        double stallAoACn;
        double stallAoANCn;
        double cdminAoA;
        double cdmin;
        std::vector<double> aoa;
        std::vector<double> cl;
        std::vector<double> cd;
        std::vector<double> cm;
    };

    // 读取翼型数据
    // 读取翼型数据
    // listFilename defaults to empty, in which case it looks for "airfoil_list.txt" in dataDir? 
    // No, let's make it explicit or use a distinct overload. 
    // Simplest: readAirfoils(dataDir, listFile)
    std::vector<AirfoilData> readAirfoils(
        const std::string &dataDir,
        const std::string &listFilename, 
        bool verbose = false);

} // namespace fvw

#endif // FVW_AIRFOIL_H