#pragma once

#include <string>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <optional>

// Simple wrapper for command line configuration
class ConfigManager {
private:
    nlohmann::json m_config;
    bool m_isLoaded = false;

public:
    ConfigManager() = default;

    // Load configuration from file
    bool loadFromFile(const std::string& filename) {
        try {
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Error: Unable to open config file: " << filename << std::endl;
                return false;
            }

            file >> m_config;
            m_isLoaded = true;
            return true;
        }
        catch (const std::exception& e) {
            std::cerr << "Error parsing config file: " << e.what() << std::endl;
            return false;
        }
    }

    // Save configuration to file
    bool saveToFile(const std::string& filename) const {
        try {
            std::ofstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
                return false;
            }

            file << m_config.dump(4); // Pretty print with 4 spaces
            return true;
        }
        catch (const std::exception& e) {
            std::cerr << "Error saving config file: " << e.what() << std::endl;
            return false;
        }
    }

    // Get entire config as JSON
    const nlohmann::json& getConfig() const {
        return m_config;
    }

    // Check if a key exists
    bool hasKey(const std::string& key) const {
        return m_config.contains(key);
    }

    // Get a value with type checking and default value
    template<typename T>
    T getValue(const std::string& key, const T& defaultValue) const {
        if (!m_config.contains(key)) {
            return defaultValue;
        }

        try {
            return m_config[key].get<T>();
        }
        catch (const std::exception& e) {
            std::cerr << "Error getting value for key '" << key << "': " << e.what() << std::endl;
            return defaultValue;
        }
    }

    // Set a value
    template<typename T>
    void setValue(const std::string& key, const T& value) {
        m_config[key] = value;
    }

    // Get a nested value using dot notation (e.g., "ga.population_size")
    template<typename T>
    std::optional<T> getNestedValue(const std::string& path) const {
        std::vector<std::string> parts;

        // Split path by dots
        size_t start = 0;
        size_t end = path.find('.');
        while (end != std::string::npos) {
            parts.push_back(path.substr(start, end - start));
            start = end + 1;
            end = path.find('.', start);
        }
        parts.push_back(path.substr(start));

        // Navigate through the JSON
        nlohmann::json current = m_config;
        for (size_t i = 0; i < parts.size() - 1; ++i) {
            if (!current.contains(parts[i])) {
                return std::nullopt;
            }
            current = current[parts[i]];
        }

        // Get the final value
        if (current.contains(parts.back())) {
            try {
                return current[parts.back()].get<T>();
            }
            catch (const std::exception& e) {
                std::cerr << "Error getting nested value for path '" << path << "': " << e.what() << std::endl;
                return std::nullopt;
            }
        }

        return std::nullopt;
    }

    // Check if configuration is loaded
    bool isLoaded() const {
        return m_isLoaded;
    }
};