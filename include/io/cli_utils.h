#ifndef FVW_CLI_UTILS_H
#define FVW_CLI_UTILS_H

#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>

namespace fvw {
    namespace cli {
        
        // --- Colors ---
        const std::string RESET   = "\033[0m";
        const std::string RED     = "\033[31m";
        const std::string GREEN   = "\033[32m";
        const std::string YELLOW  = "\033[33m";
        const std::string BLUE    = "\033[34m";
        const std::string MAGENTA = "\033[35m";
        const std::string CYAN    = "\033[36m";
        const std::string WHITE   = "\033[37m";
        const std::string BOLD    = "\033[1m";

        inline void print_header(const std::string& title) {
            std::cout << "\n" << BOLD << CYAN << "=== " << title << " ===" << RESET << "\n";
        }

        inline void print_info(const std::string& key, const std::string& value) {
            std::cout << GREEN << "  [+]" << RESET << " " << std::setw(20) << std::left << key << ": " << value << "\n";
        }

        inline void print_warning(const std::string& msg) {
            std::cout << YELLOW << "[WARNING] " << msg << RESET << "\n";
        }

        inline void print_error(const std::string& msg) {
            std::cerr << RED << "[ERROR] " << msg << RESET << "\n";
        }

        // --- Progress Bar ---
        class ProgressBar {
        public:
            ProgressBar(int total, int width = 50) : m_total(total), m_width(width) {}

            void update(int current, double elapsed_sec = 0.0) {
                if (m_total <= 0) return;
                
                float progress = static_cast<float>(current) / m_total;
                int pos = static_cast<int>(m_width * progress);

                // Carraige return to overwrite line
                std::cout << "\r" << "[";
                for (int i = 0; i < m_width; ++i) {
                    if (i < pos) std::cout << "=";
                    else if (i == pos) std::cout << ">";
                    else std::cout << " ";
                }
                std::cout << "] " << int(progress * 100.0) << "% ";
                
                // Info
                if (elapsed_sec > 0) {
                    std::cout << std::fixed << std::setprecision(1) << "Time: " << elapsed_sec << "s";
                }
                
                std::cout << std::flush;
            }

            void finish() {
                std::cout << "\n" << RESET;
            }

        private:
            int m_total;
            int m_width;
        };

    } // namespace cli
} // namespace fvw

#endif // FVW_CLI_UTILS_H
