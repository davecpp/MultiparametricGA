# Multi-Parametric VLSI Placement Optimization

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A C++ framework for optimizing VLSI (Very Large Scale Integration) chip layout placement using genetic algorithms with multi-parametric optimization and distributed computation.

## Overview

This project implements an advanced multi-objective optimization system for integrated circuit (IC) placement. The framework uses evolutionary algorithms to find optimal placement configurations that balance various competing design parameters:

- **Wire length optimization** - minimizing the total connection length between cells
- **Thermal distribution** - ensuring even heat distribution across the chip
- **Design Rule Checking (DRC)** - maintaining minimum distance requirements and avoiding overlaps
- **Power consumption** - reducing both dynamic and static power consumption
- **Parasitic effects** - minimizing cross-talk and signal integrity issues

## Architecture

The system consists of several integrated subsystems:

### Core Components

- **Genetic Algorithm Engine** - Handles population evolution with crossover, mutation, and selection operators
- **Multi-Parametric Fitness Function** - Combines multiple objectives with weighted coefficients
- **Geometric Operations Module** - Manages geometric calculations and spatial queries
- **Design Rule Checker** - Validates layout against design constraints
- **Penalty System** - Efficiently identifies and penalizes rule violations
- **Specialized Sub-Populations** - Optimizes different aspects of the solution independently
- **Distributed Computation System** - Parallelizes calculations across multiple threads or nodes

### Key Files

- `GeneticAlgorithm.h` - Main evolutionary algorithm implementation
- `MultiParametricFitness.h` - Multi-objective fitness function
- `DRCChecker.h` - Design rule checking system
- `Chromosome.h` - Represents placement solutions
- `Cell.h` - Represents IC components
- `Scheme.h` - Defines the overall circuit scheme
- `GeoJSONExporter.h` - Exports results to GeoJSON format
- `Population.h` - Manages collections of chromosomes

## Features

### Advanced Genetic Algorithm

- Adaptive crossover and mutation rates
- Tournament and roulette wheel selection
- Elitism to preserve best solutions
- Dynamic convergence criteria

### Specialized Sub-Populations

The framework uses parallel sub-populations, each focused on optimizing different criteria:

- **WireLength** - Wire length only
- **Thermal** - Thermal distribution only
- **DRC** - Design Rule Checking only
- **Power** - Power consumption only
- **Parasitic** - Parasitic effects only
- **WireThermal** - Wire length + thermal
- **PowerDRC** - Power + DRC
- **Full** - All parameters

### Efficient DRC Penalty System

Rather than integrating DRC into the fitness function (which would be computation-intensive), the system uses an efficient post-filtering approach:

- Evaluates DRC violations after each generation
- Ranks solutions by violation severity
- Penalizes or eliminates worst offenders
- Adapts penalty coefficients dynamically

### Distributed Computation

- Parallel evolution of sub-populations
- Multi-threaded fitness evaluation
- Optional cluster support for large-scale optimization

## Usage

### Basic Example

```cpp
// Create and configure genetic algorithm
GeneticAlgorithm ga;
GeneticAlgorithm::Parameters params;
params.populationSize = 100;
params.maxIterations = 200;
params.useSpecializedPopulations = true;
ga.setParameters(params);

// Configure fitness function weights
MultiParametricFitness::Weights weights;
weights.wireLength = 1.0;
weights.thermal = 0.5;
weights.drc = 2.0;
weights.power = 0.3;
weights.parasitic = 0.2;

// Configure DRC parameters
DRCChecker::DRCParams drcParams;
drcParams.minDistance = 2;
drcParams.overlapPenalty = 10.0;
drcParams.distancePenalty = 5.0;

// Set up fitness and DRC
ga.setFitnessParameters(weights, thermalParams, powerParams);
ga.setDRCParameters(drcParams);

// Load scheme
Scheme scheme;
scheme.loadFromFile("my_circuit.json");

// Run optimization
Chromosome bestSolution = ga.run(scheme);

// Export results
GeoJSONExporter::exportToFile(bestSolution, scheme, "placement_result.json");
```

### Command Line Interface

The system also provides a command-line interface for easier integration into design flows:

```bash
./vlsi_optimizer --config=config.json --input=circuit.json --output=result.json
```

## Configuration

The genetic algorithm can be configured via JSON:

```json
{
  "population_size": 100,
  "max_iterations": 100,
  "convergence_threshold": 0.001,
  "filter_coefficient": 0.4,
  "mutation_probability": 0.1,
  "drc_filter_percentage": 0.1,
  "drc_check_interval": 10,
  "use_specialized_populations": true,
  "specialized_evolutions": 20,
  "specialized_best_percentage": 0.2,
  "specialized_random_percentage": 0.05
}
```

## Building

### Requirements

- C++17 compatible compiler
- CMake 3.10+
- nlohmann/json library

### Build Steps

```bash
mkdir build
cd build
cmake ..
make
```

## Documentation

For more detailed information about the algorithms and implementation details, refer to the documentation in the `docs` directory.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The project draws inspiration from research on multi-parametric optimization in VLSI design
- Special thanks to contributors and reviewers
