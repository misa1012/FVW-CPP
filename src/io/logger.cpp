#include "io/logger.h"
#include "io/cli_utils.h"
#include <iomanip>
#include <chrono>
#include <ctime>

namespace fvw {

    std::ofstream Logger::m_logFile;
    std::mutex Logger::m_mutex;

    void Logger::init(const std::string& filepath) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_logFile.is_open()) {
            m_logFile.close();
        }
        m_logFile.open(filepath);
        if (m_logFile.is_open()) {
            // Write header
            auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            m_logFile << "=== FVW-CPP Simulation Log - Started at " << std::ctime(&now); // ctime adds newline
            m_logFile << "========================================================\n" << std::endl;
        } else {
            std::cerr << "Failed to open log file: " << filepath << std::endl;
        }
    }

    void Logger::close() {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_logFile.is_open()) {
            m_logFile << "\n=== Log Closed ===" << std::endl;
            m_logFile.close();
        }
    }

    void Logger::log(const std::string& message) {
        log("LOG", message); // Default tag handling
    }

    void Logger::log(const std::string& tag, const std::string& message) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_logFile.is_open()) {
            // Format: "  |-- [TAG]    Msg"
            // Pad tag to fixed width (e.g., 8 chars)
            std::string paddedTag = "[" + tag + "]";
            while(paddedTag.length() < 10) paddedTag += " ";
            
            m_logFile << "  |-- " << paddedTag << message << std::endl;
        }
    }

    void Logger::section_header(double time, int step) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_logFile.is_open()) {
            m_logFile << "\nTime = " << std::fixed << std::setprecision(3) << time << " s (Step " << step << ")" << std::endl;
            m_logFile << "--------------------------------------------------------" << std::endl;
        }
    }

    void Logger::info(const std::string& key, const std::string& value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        // Console: Pretty print
        cli::print_info(key, value);
        // File: Simple print
        if (m_logFile.is_open()) {
            m_logFile << "[INFO] " << key << ": " << value << std::endl;
        }
    }

    void Logger::info(const std::string& message) {
        std::lock_guard<std::mutex> lock(m_mutex);
        // Console
        std::cout << message << std::endl;
        // File
        if (m_logFile.is_open()) {
            m_logFile << "[INFO] " << message << std::endl;
        }
    }

    void Logger::warn(const std::string& message) {
        std::lock_guard<std::mutex> lock(m_mutex);
        // Console: Yellow
        cli::print_warning(message);
        // File
        if (m_logFile.is_open()) {
            m_logFile << "[WARN] " << message << std::endl;
        }
    }

    void Logger::error(const std::string& message) {
        std::lock_guard<std::mutex> lock(m_mutex);
        // Console: Red
        cli::print_error(message);
        // File
        if (m_logFile.is_open()) {
            m_logFile << "[ERR ] " << message << std::endl;
        }
    }

} // namespace fvw
