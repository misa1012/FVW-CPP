#include "airfoil.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace fvw
{

    std::vector<AirfoilData> readAirfoils(const std::string &dataDir, const std::string &listFilename, bool verbose)
    {
        std::vector<std::string> airfoilNames;
        
        // Read names from list file
        std::string listPath = listFilename;
        // If listFilename is relative and doesn't exist, maybe check if it's relative to dataDir?
        // But for simplicity, we expect listFilename to be a valid path (absolute or relative to cwd).
        // Let's assume listFilename is the full path to the list file.
        
        std::ifstream listFile(listPath);
        if (!listFile.is_open()) {
             // Try combining dataDir + listFilename if simple open fails? 
             // Let's keep it simple: user provides path. 
             std::cerr << "Error: Could not open airfoil list file: " << listPath << std::endl;
             // Fallback to hardcoded for now? No, we want to remove hardcoding.
             return {};
        }
        
        std::string name;
        while (std::getline(listFile, name)) {
            // Trim whitespace
            name.erase(0, name.find_first_not_of(" \t\r\n"));
            name.erase(name.find_last_not_of(" \t\r\n") + 1);
            if (!name.empty()) {
                airfoilNames.push_back(name);
            }
        }
        listFile.close();

        std::vector<AirfoilData> airfoils(airfoilNames.size());
        for (size_t i = 0; i < airfoilNames.size(); ++i)
        {
            // Ensure dataDir has trailing slash if not empty
            std::string dir = dataDir;
            if (!dir.empty() && dir.back() != '/') dir += "/";
            
            std::string filePath = dir + airfoilNames[i] + ".dat";
            if (verbose)
            {
                std::cout << "Reading: " << filePath << std::endl;
            }
            std::ifstream file(filePath);
            if (!file.is_open())
            {
                std::cerr << "Error opening " << filePath << std::endl;
                continue;
            }

            std::vector<std::string> header(14);
            std::string line;
            for (int j = 0; j < 14; ++j)
            {
                if (!std::getline(file, line))
                {
                    std::cerr << "Error: " << filePath << " has fewer than 14 header lines" << std::endl;
                    file.close();
                    break;
                }
                header[j] = line;
            }

            try
            {
                std::istringstream iss4(header[4]);
                iss4 >> airfoils[i].stallAoA;
                std::istringstream iss8(header[8]);
                iss8 >> airfoils[i].cn0AoA;
                std::istringstream iss9(header[9]);
                iss9 >> airfoils[i].lift0Cn;
                std::istringstream iss10(header[10]);
                iss10 >> airfoils[i].stallAoACn;
                std::istringstream iss11(header[11]);
                iss11 >> airfoils[i].stallAoANCn;
                std::istringstream iss12(header[12]);
                iss12 >> airfoils[i].cdminAoA;
                std::istringstream iss13(header[13]);
                iss13 >> airfoils[i].cdmin;
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error parsing header in " << filePath << ": " << e.what() << std::endl;
                file.close();
                continue;
            }

            while (std::getline(file, line))
            {
                line.erase(0, line.find_first_not_of(" \t\r\n"));
                line.erase(line.find_last_not_of(" \t\r\n") + 1);
                if (line.empty())
                    continue;
                std::istringstream iss(line);
                double aoa, cl, cd, cm;
                if (iss >> aoa >> cl >> cd >> cm)
                {
                    airfoils[i].aoa.push_back(aoa);
                    airfoils[i].cl.push_back(cl);
                    airfoils[i].cd.push_back(cd);
                    airfoils[i].cm.push_back(cm);
                }
                else
                {
                    std::cerr << "Warning: Invalid data line in " << filePath << ": '" << line << "'" << std::endl;
                }
            }
            file.close();
            if (verbose)
            {
                std::cout << "Parsed " << filePath << ", points: " << airfoils[i].aoa.size() << std::endl;
            }
        }
        return airfoils;
    }

} // namespace fvw