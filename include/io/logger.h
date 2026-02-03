#ifndef FVW_LOGGER_H
#define FVW_LOGGER_H

#include <string>
#include <fstream>
#include <iostream>
#include <mutex>

namespace fvw {

    class Logger {
    public:
        // Initialize the logger with a file path
        static void init(const std::string& filepath);

        // Close the log file
        static void close();

        // Log with a specific TAG (e.g., [WAKE], [IO]) and auto-indentation
        // Format: "  |-- [TAG]    Message"
        static void log(const std::string& tag, const std::string& message);
        static void log(const std::string& message); // Legacy overload for simple messages

        // Log a section header for a new timestep
        // Format: "Time = 0.50 s (Step 10)"
        static void section_header(double time, int step);

        // Log to FILE and CONSOLE (Standard info, e.g., "Case started")
        static void info(const std::string& key, const std::string& value);
        static void info(const std::string& message);

        // Log to FILE and CONSOLE (Warnings, e.g., "Convergence failure")
        static void warn(const std::string& message);

        // Log to FILE and CONSOLE (Errors)
        static void error(const std::string& message);

    private:
        static std::ofstream m_logFile;
        static std::mutex m_mutex;
    };

} // namespace fvw

#endif // FVW_LOGGER_H
