#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <optional>

// Class for parsing command line arguments
class CommandLineParser {
private:
    std::unordered_map<std::string, std::string> m_args;
    std::vector<std::string> m_positionalArgs;

public:
    // Parse command line arguments
    void parse(int argc, char* argv[]) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            // Check if it's a named argument
            if (arg.substr(0, 2) == "--") {
                // Format: --key=value or --key value
                size_t equalPos = arg.find('=');
                if (equalPos != std::string::npos) {
                    // --key=value format
                    std::string key = arg.substr(2, equalPos - 2);
                    std::string value = arg.substr(equalPos + 1);
                    m_args[key] = value;
                }
                else {
                    // --key value format (if next arg exists and doesn't start with --)
                    std::string key = arg.substr(2);
                    if (i + 1 < argc && argv[i + 1][0] != '-') {
                        m_args[key] = argv[i + 1];
                        ++i; // Skip the value in the next iteration
                    }
                    else {
                        // Boolean flag (--flag with no value)
                        m_args[key] = "true";
                    }
                }
            }
            else if (arg[0] == '-') {
                // Short option (-k value or -k=value)
                if (arg.length() > 1) {
                    if (arg.length() > 2 && arg[2] == '=') {
                        // -k=value format
                        std::string key = arg.substr(1, 1);
                        std::string value = arg.substr(3);
                        m_args[key] = value;
                    }
                    else {
                        // -k value format
                        std::string key = arg.substr(1, 1);
                        if (i + 1 < argc && argv[i + 1][0] != '-') {
                            m_args[key] = argv[i + 1];
                            ++i; // Skip the value in the next iteration
                        }
                        else {
                            // Boolean flag (-f with no value)
                            m_args[key] = "true";
                        }
                    }
                }
            }
            else {
                // Positional argument
                m_positionalArgs.push_back(arg);
            }
        }
    }

    // Get a named argument with a default value
    std::string getArg(const std::string& key, const std::string& defaultValue = "") const {
        auto it = m_args.find(key);
        return (it != m_args.end()) ? it->second : defaultValue;
    }

    // Get an optional named argument
    std::optional<std::string> getOptionalArg(const std::string& key) const {
        auto it = m_args.find(key);
        return (it != m_args.end()) ? std::optional<std::string>(it->second) : std::nullopt;
    }

    // Get a positional argument
    std::string getPositionalArg(size_t index, const std::string& defaultValue = "") const {
        return (index < m_positionalArgs.size()) ? m_positionalArgs[index] : defaultValue;
    }

    // Get all positional arguments
    const std::vector<std::string>& getPositionalArgs() const {
        return m_positionalArgs;
    }

    // Check if a named argument exists
    bool hasArg(const std::string& key) const {
        return m_args.find(key) != m_args.end();
    }

    // Print help message for the application
    static void printHelp(const std::string& programName) {
        std::cout << "Usage: " << programName << " [OPTIONS]\n\n"
            << "Options:\n"
            << "  --config=<file>          Path to configuration file\n"
            << "  --input=<file>           Path to input scheme file\n"
            << "  --output=<file>          Path to output file\n"
            << "  --help                   Show this help message\n"
            << std::endl;
    }
};