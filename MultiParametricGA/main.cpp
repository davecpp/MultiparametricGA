#include "Cell.h"
#include "Chromosome.h"
#include "Coordinate.h"
#include "GeneticAlgorithm.h"
#include "MultiParametricFitness.h"
#include "DRCChecker.h" // Добавляем заголовок DRCChecker
#include "Scheme.h"
#include "ConfigParser.h"
#include "CommandLine.h"
#include "GeoJSONExporter.h"
#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <ctime>

// Default configuration file
constexpr const char* DEFAULT_CONFIG_FILE = "config.json";

// Create a default configuration
nlohmann::json createDefaultConfig() {
    nlohmann::json config;

    // GA parameters
    config["ga"] = {
        {"population_size", 100},
        {"max_iterations", 200},
        {"convergence_threshold", 0.001},
        {"filter_coefficient", 0.4},
        {"mutation_probability", 0.1},
        {"drc_filter_percentage", 0.1},
        {"drc_check_interval", 10},
        {"use_specialized_populations", true},
        {"specialized_evolutions", 20},
        {"specialized_best_percentage", 0.2},
        {"specialized_random_percentage", 0.05}
    };

    // Fitness weights
    config["fitness"] = {
        {"weights", {
            {"wire_length", 1.0},
            {"thermal", 0.5},
            {"drc", 2.0},
            {"power", 0.3},
            {"parasitic", 0.2}
        }},
        {"thermal", {
            {"max_allowed_temp", 100.0},
            {"ambient_temp", 25.0},
            {"thermal_conductivity", 150.0},
            {"power_density", 0.1}
        }},
        {"power", {
            {"cell_power", 0.1},
            {"wire_capacitance", 0.02},
            {"voltage_factor", 1.0},
            {"frequency", 1.0},
            {"leakage_factor", 0.05}
        }}
    };

    // DRC parameters
    config["drc"] = {
        {"min_distance", 2},
        {"overlap_penalty", 10.0},
        {"distance_penalty", 5.0}
    };

    return config;
}

// Load configuration from file or create default
ConfigManager loadConfiguration(const std::string& configFile) {
    ConfigManager config;

    // Try to load from file
    if (!configFile.empty() && config.loadFromFile(configFile)) {
        std::cout << "Configuration loaded from " << configFile << std::endl;
    }
    else {
        // Create default configuration
        std::cout << "Using default configuration" << std::endl;
        nlohmann::json defaultConfig = createDefaultConfig();
        config = ConfigManager();
        config.setValue("ga", defaultConfig["ga"]);
        config.setValue("fitness", defaultConfig["fitness"]);
        config.setValue("drc", defaultConfig["drc"]);

        // Save default configuration if file was specified
        if (!configFile.empty()) {
            config.saveToFile(configFile);
            std::cout << "Default configuration saved to " << configFile << std::endl;
        }
    }

    return config;
}

// Configure genetic algorithm from config
void configureGA(GeneticAlgorithm& ga, const ConfigManager& config) {
    // Set GA parameters
    GeneticAlgorithm::Parameters params;
    const auto& gaConfig = config.getConfig()["ga"];

    params.populationSize = gaConfig.value("population_size", 100);
    params.maxIterations = gaConfig.value("max_iterations", 200);
    params.convergenceThreshold = gaConfig.value("convergence_threshold", 0.001);
    params.filterCoefficient = gaConfig.value("filter_coefficient", 0.4);
    params.mutationProbability = gaConfig.value("mutation_probability", 0.1);
    params.drcFilterPercentage = gaConfig.value("drc_filter_percentage", 0.1);
    params.drcCheckInterval = gaConfig.value("drc_check_interval", 10);
    params.useSpecializedPopulations = gaConfig.value("use_specialized_populations", true);
    params.specializedEvolutions = gaConfig.value("specialized_evolutions", 20);
    params.specializedBestPercentage = gaConfig.value("specialized_best_percentage", 0.2);
    params.specializedRandomPercentage = gaConfig.value("specialized_random_percentage", 0.05);

    ga.setParameters(params);

    // Set fitness parameters
    const auto& fitnessConfig = config.getConfig()["fitness"];

    // Weights
    MultiParametricFitness::Weights weights;
    const auto& weightsConfig = fitnessConfig["weights"];
    weights.wireLength = weightsConfig.value("wire_length", 1.0);
    weights.thermal = weightsConfig.value("thermal", 0.5);
    weights.drc = weightsConfig.value("drc", 2.0);
    weights.power = weightsConfig.value("power", 0.3);
    weights.parasitic = weightsConfig.value("parasitic", 0.2);

    // Thermal parameters
    MultiParametricFitness::ThermalModelParams thermalParams;
    const auto& thermalConfig = fitnessConfig["thermal"];
    thermalParams.maxAllowedTemp = thermalConfig.value("max_allowed_temp", 100.0);
    thermalParams.ambientTemp = thermalConfig.value("ambient_temp", 25.0);
    thermalParams.thermalConductivity = thermalConfig.value("thermal_conductivity", 150.0);
    thermalParams.powerDensity = thermalConfig.value("power_density", 0.1);

    // Power parameters
    MultiParametricFitness::PowerModelParams powerParams;
    const auto& powerConfig = fitnessConfig["power"];
    powerParams.cellPower = powerConfig.value("cell_power", 0.1);
    powerParams.wireCapacitance = powerConfig.value("wire_capacitance", 0.02);
    powerParams.voltageFactor = powerConfig.value("voltage_factor", 1.0);
    powerParams.frequency = powerConfig.value("frequency", 1.0);
    powerParams.leakageFactor = powerConfig.value("leakage_factor", 0.05);

    ga.setFitnessParameters(weights, thermalParams, powerParams);

    // Set DRC parameters
    const auto& drcConfig = config.getConfig()["drc"];
    DRCChecker::DRCParams drcParams;
    drcParams.minDistance = drcConfig.value("min_distance", 2);
    drcParams.overlapPenalty = drcConfig.value("overlap_penalty", 10.0);
    drcParams.distancePenalty = drcConfig.value("distance_penalty", 5.0);

    ga.setDRCParameters(drcParams);
}

// Generate test scheme for testing purposes
Scheme generateTestScheme(int cellCount, int fieldRows, int fieldCols) {
    return Scheme::generateTestScheme(cellCount, fieldRows, fieldCols);
}

// Main function
int main(int argc, char* argv[]) {
    // Initialize random seed
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // Parse command line arguments
    CommandLineParser cmdParser;
    cmdParser.parse(argc, argv);

    // Check for help flag
    if (cmdParser.hasArg("help") || cmdParser.hasArg("h")) {
        CommandLineParser::printHelp(argv[0]);
        return 0;
    }

    // Get file paths
    std::string configFile = cmdParser.getArg("config", DEFAULT_CONFIG_FILE);
    std::string inputFile = cmdParser.getArg("input", "");
    std::string outputFile = cmdParser.getArg("output", "placement_output.geojson");

    // Load configuration
    ConfigManager config = loadConfiguration(configFile);

    // Create genetic algorithm
    GeneticAlgorithm ga;
    configureGA(ga, config);

    // Set verbosity
    ga.setVerbose(cmdParser.hasArg("verbose") || cmdParser.hasArg("v"));

    // Load or generate scheme
    Scheme scheme;
    bool schemeLoaded = false;

    if (!inputFile.empty()) {
        std::cout << "Loading scheme from " << inputFile << "..." << std::endl;
        schemeLoaded = scheme.loadFromFile(inputFile);
    }

    if (!schemeLoaded) {
        // Generate a test scheme if input file not specified or loading failed
        int cellCount = cmdParser.hasArg("cells") ?
            std::stoi(cmdParser.getArg("cells")) : 100;
        int fieldRows = cmdParser.hasArg("rows") ?
            std::stoi(cmdParser.getArg("rows")) : 20;
        int fieldCols = cmdParser.hasArg("cols") ?
            std::stoi(cmdParser.getArg("cols")) : 20;

        std::cout << "Generating test scheme with " << cellCount << " cells, "
            << fieldRows << "x" << fieldCols << " grid..." << std::endl;

        scheme = generateTestScheme(cellCount, fieldRows, fieldCols);

        // Save generated scheme if requested
        if (cmdParser.hasArg("save-scheme")) {
            std::string schemeFile = cmdParser.getArg("save-scheme");
            if (scheme.saveToFile(schemeFile)) {
                std::cout << "Test scheme saved to " << schemeFile << std::endl;
            }
        }
    }

    // Check if scheme is valid
    if (!scheme.isValid()) {
        std::cerr << "Error: Invalid scheme" << std::endl;
        return 1;
    }

    // Run the genetic algorithm
    std::cout << "Running genetic algorithm..." << std::endl;
    auto startTime = std::chrono::high_resolution_clock::now();

    Chromosome bestPlacement = ga.run(scheme);

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // Print statistics
    ga.printStatistics();

    // Export placement to GeoJSON
    std::cout << "Exporting placement to " << outputFile << "..." << std::endl;
    if (GeoJSONExporter::exportToFile(bestPlacement, scheme, outputFile)) {
        std::cout << "Placement exported successfully" << std::endl;
    }
    else {
        std::cerr << "Error: Failed to export placement" << std::endl;
        return 1;
    }

    return 0;
}
