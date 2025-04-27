#pragma once

#include "Chromosome.h"
#include "Population.h"
#include "DRCChecker.h"          // Сначала включаем DRCChecker
#include "FitnessFunction.h"      // Затем FitnessFunction, который включает MultiParametricFitness
#include "Scheme.h"
#include <vector>
#include <random>
#include <thread>
#include <functional>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <memory>
#include <nlohmann/json.hpp>

// Type for specialized population type
enum class SpecializedPopType {
	WireLength,    // Wire length only
	Thermal,       // Thermal distribution only
	DRC,           // Design Rule Checking only
	Power,         // Power consumption only
	Parasitic,     // Parasitic effects only
	WireThermal,   // Wire length + thermal
	PowerDRC,      // Power + DRC
	Full           // All parameters
};

// Forward declarations
class SpecializedPopulation;
class SpecializedPopulationsManager;

// Main Genetic Algorithm class
class GeneticAlgorithm {
public:
	// Parameters for the genetic algorithm
	struct Parameters {
		size_t populationSize = 100;
		size_t maxIterations = 100;
		double convergenceThreshold = 0.001;
		double filterCoefficient = 0.4;
		double mutationProbability = 0.1;
		double drcFilterPercentage = 0.1;
		size_t drcCheckInterval = 10;
		bool useSpecializedPopulations = true;
		size_t specializedEvolutions = 20;
		double specializedBestPercentage = 0.2;
		double specializedRandomPercentage = 0.05;

		// Convert to JSON
		nlohmann::json toJson() const {
			nlohmann::json j;
			j["population_size"] = populationSize;
			j["max_iterations"] = maxIterations;
			j["convergence_threshold"] = convergenceThreshold;
			j["filter_coefficient"] = filterCoefficient;
			j["mutation_probability"] = mutationProbability;
			j["drc_filter_percentage"] = drcFilterPercentage;
			j["drc_check_interval"] = drcCheckInterval;
			j["use_specialized_populations"] = useSpecializedPopulations;
			j["specialized_evolutions"] = specializedEvolutions;
			j["specialized_best_percentage"] = specializedBestPercentage;
			j["specialized_random_percentage"] = specializedRandomPercentage;
			return j;
		}

		// Load from JSON
		static Parameters fromJson(const nlohmann::json& j) {
			Parameters params;

			if (j.contains("population_size")) {
				params.populationSize = j["population_size"].get<size_t>();
			}

			if (j.contains("max_iterations")) {
				params.maxIterations = j["max_iterations"].get<size_t>();
			}

			if (j.contains("convergence_threshold")) {
				params.convergenceThreshold = j["convergence_threshold"].get<double>();
			}

			if (j.contains("filter_coefficient")) {
				params.filterCoefficient = j["filter_coefficient"].get<double>();
			}

			if (j.contains("mutation_probability")) {
				params.mutationProbability = j["mutation_probability"].get<double>();
			}

			if (j.contains("drc_filter_percentage")) {
				params.drcFilterPercentage = j["drc_filter_percentage"].get<double>();
			}

			if (j.contains("drc_check_interval")) {
				params.drcCheckInterval = j["drc_check_interval"].get<size_t>();
			}

			if (j.contains("use_specialized_populations")) {
				params.useSpecializedPopulations = j["use_specialized_populations"].get<bool>();
			}

			if (j.contains("specialized_evolutions")) {
				params.specializedEvolutions = j["specialized_evolutions"].get<size_t>();
			}

			if (j.contains("specialized_best_percentage")) {
				params.specializedBestPercentage = j["specialized_best_percentage"].get<double>();
			}

			if (j.contains("specialized_random_percentage")) {
				params.specializedRandomPercentage = j["specialized_random_percentage"].get<double>();
			}

			return params;
		}
	};

	// Statistics about the algorithm run
	struct Statistics {
		double initialFitness = 0.0;
		double finalFitness = 0.0;
		double bestWireLengthFitness = 0.0;
		double bestThermalFitness = 0.0;
		double bestDRCFitness = 0.0;
		double bestPowerFitness = 0.0;
		double bestParasiticFitness = 0.0;
		size_t iterations = 0;
		std::chrono::milliseconds executionTime{ 0 };

		// Convert to JSON
		nlohmann::json toJson() const {
			nlohmann::json j;
			j["initial_fitness"] = initialFitness;
			j["final_fitness"] = finalFitness;
			j["best_wire_length_fitness"] = bestWireLengthFitness;
			j["best_thermal_fitness"] = bestThermalFitness;
			j["best_drc_fitness"] = bestDRCFitness;
			j["best_power_fitness"] = bestPowerFitness;
			j["best_parasitic_fitness"] = bestParasiticFitness;
			j["iterations"] = iterations;
			j["execution_time_ms"] = executionTime.count();
			return j;
		}
	};

private:
	Parameters m_params;
	MultiParametricFitness m_fitness;
	DRCChecker m_drcChecker;
	Statistics m_stats;
	bool m_verbose = false;

public:
	// Constructor
	GeneticAlgorithm() = default;

	GeneticAlgorithm(const Parameters& params)
		: m_params(params) {
	}

	// Set parameters
	void setParameters(const Parameters& params) {
		m_params = params;
	}

	// Set fitness parameters
	void setFitnessParameters(const MultiParametricFitness::Weights& weights,
		const MultiParametricFitness::ThermalModelParams& thermalParams,
		const MultiParametricFitness::PowerModelParams& powerParams) {
		m_fitness.setWeights(weights);
		m_fitness.setThermalParams(thermalParams);
		m_fitness.setPowerParams(powerParams);
		m_fitness.setDRCChecker(&m_drcChecker);
	}

	// Set DRC parameters
	void setDRCParameters(const DRCChecker::DRCParams& params) {
		m_drcChecker.setParams(params);
	}

	// Set verbosity
	void setVerbose(bool verbose) {
		m_verbose = verbose;
	}

	// Get statistics
	const Statistics& getStatistics() const {
		return m_stats;
	}

	// Print statistics
	void printStatistics() const {
		std::cout << "\n------ Genetic Algorithm Statistics ------\n";
		std::cout << std::fixed << std::setprecision(4);
		std::cout << "Initial fitness: " << m_stats.initialFitness << "\n";
		std::cout << "Final fitness: " << m_stats.finalFitness << "\n";
		std::cout << "Improvement: " <<
			(m_stats.initialFitness > 0 ?
				(1.0 - m_stats.finalFitness / m_stats.initialFitness) * 100.0 : 0.0) <<
			"%\n\n";

		std::cout << "Best solution components:\n";
		std::cout << "  Wire length fitness: " << m_stats.bestWireLengthFitness << "\n";
		std::cout << "  Thermal fitness: " << m_stats.bestThermalFitness << "\n";
		std::cout << "  DRC fitness: " << m_stats.bestDRCFitness << "\n";
		std::cout << "  Power fitness: " << m_stats.bestPowerFitness << "\n";
		std::cout << "  Parasitic fitness: " << m_stats.bestParasiticFitness << "\n\n";

		std::cout << "Iterations: " << m_stats.iterations << "\n";
		std::cout << "Execution time: " << m_stats.executionTime.count() << " ms\n";
		std::cout << "-----------------------------------------\n\n";
	}

	// Main run method
	Chromosome run(const Scheme& scheme);

	// Run the genetic algorithm with specialized populations
	Chromosome runWithSpecializedPopulations(const Scheme& scheme);

	// Run standard genetic algorithm
	Chromosome runStandardEvolution(Population& population, const Scheme& scheme);

	// Static utility methods

	// Crossover function
	static Chromosome crossover(const Chromosome& parent1,
		const Chromosome& parent2,
		const Scheme& scheme);

	// Mutation function
	static void mutate(Chromosome& chromosome, const Scheme& scheme);
};

// Class for specialized populations
class SpecializedPopulation {
private:
	Population m_population;
	MultiParametricFitness m_fitness;
	SpecializedPopType m_type;
	size_t m_evolutions = 20;

public:
	// Constructor
	SpecializedPopulation(size_t size, SpecializedPopType type)
		: m_population(size), m_type(type) {
		// Set up fitness function based on type
		configureForType(type);
	}

	// Configure fitness function for specific type
	void configureForType(SpecializedPopType type);

	// Set evolution iterations
	void setEvolutions(size_t evolutions) {
		m_evolutions = evolutions;
	}

	// Set fitness function parameters
	void setFitnessParameters(const MultiParametricFitness::ThermalModelParams& thermalParams,
		const MultiParametricFitness::PowerModelParams& powerParams,
		DRCChecker* drcChecker);

	// Initialize the population
	void initialize(const Scheme& scheme);

	// Evolve the population
	void evolve(const Scheme& scheme,
		const std::function<Chromosome(const Chromosome&, const Chromosome&, const Scheme&)>& crossoverFunc,
		const std::function<void(Chromosome&, const Scheme&)>& mutateFunc,
		double mutationProbability);

	// Get best individuals
	std::vector<Individual> getBestIndividuals(double percentage) const;

	// Get random individuals
	std::vector<Individual> getRandomIndividuals(double percentage) const;
};

// Manager for specialized populations
class SpecializedPopulationsManager {
private:
	std::vector<SpecializedPopulation> m_populations;
	size_t m_basePopulationSize;

public:
	// Constructor
	SpecializedPopulationsManager(size_t basePopSize)
		: m_basePopulationSize(basePopSize) {
	}

	// Add a specialized population
	void addSpecializedPopulation(SpecializedPopType type) {
		m_populations.emplace_back(m_basePopulationSize, type);
	}

	// Set evolution iterations for all populations
	void setEvolutions(size_t evolutions) {
		for (auto& pop : m_populations) {
			pop.setEvolutions(evolutions);
		}
	}

	// Set fitness parameters for all populations
	void setFitnessParameters(const MultiParametricFitness::ThermalModelParams& thermalParams,
		const MultiParametricFitness::PowerModelParams& powerParams,
		DRCChecker* drcChecker);

	// Initialize all populations
	void initializeAll(const Scheme& scheme);

	// Evolve all populations in parallel
	void evolveAll(const Scheme& scheme,
		const std::function<Chromosome(const Chromosome&, const Chromosome&, const Scheme&)>& crossoverFunc,
		const std::function<void(Chromosome&, const Scheme&)>& mutateFunc,
		double mutationProbability);

	// Create a merged population from the best individuals
	Population createMergedPopulation(const Scheme& scheme,
		const MultiParametricFitness& fitness,
		double bestPercentage,
		double randomPercentage);

	// Get population count
	size_t getPopulationCount() const {
		return m_populations.size();
	}
};

// Implementation of main run method
inline Chromosome GeneticAlgorithm::run(const Scheme& scheme) {
	if (m_params.useSpecializedPopulations) {
		if (m_verbose) {
			std::cout << "Running with specialized populations...\n";
		}
		return runWithSpecializedPopulations(scheme);
	}
	else {
		if (m_verbose) {
			std::cout << "Running standard genetic algorithm...\n";
		}
		Population population(m_params.populationSize);
		return runStandardEvolution(population, scheme);
	}
}

// Implementation of the crossover function
inline Chromosome GeneticAlgorithm::crossover(const Chromosome& parent1,
	const Chromosome& parent2,
	const Scheme& scheme) {
	// Create a new chromosome
	Chromosome child(scheme);

	// Mark positions that are already used
	const auto& field = scheme.getField();
	Matrix<bool> usedPositions(field.rows, field.cols, false);

	// Calculate probabilities based on fitness
	// This is a placeholder - in a real implementation, use actual fitness values
	double parent1Prob = 0.5;

	// Determine gene sources based on probability
	std::vector<bool> geneFromParent1(parent1.size());
	std::random_device rd;
	std::mt19937 gen(rd());
	std::bernoulli_distribution dist(parent1Prob);

	for (size_t i = 0; i < geneFromParent1.size(); ++i) {
		geneFromParent1[i] = dist(gen);
	}

	// First, place genes from parent1 according to probability
	for (size_t i = 0; i < parent1.size(); ++i) {
		if (geneFromParent1[i]) {
			const Coordinate& pos = parent1[i];

			// Check if position is already used
			if (pos.isValid() && pos.y() < field.rows && pos.x() < field.cols &&
				!usedPositions[pos.y()][pos.x()]) {
				child[i] = pos;
				usedPositions[pos.y()][pos.x()] = true;
			}
		}
	}

	// Then, place genes from parent2 where not already placed
	for (size_t i = 0; i < parent2.size(); ++i) {
		if (!geneFromParent1[i] || !child[i].isValid()) {
			const Coordinate& pos = parent2[i];

			// Check if position is already used
			if (pos.isValid() && pos.y() < field.rows && pos.x() < field.cols &&
				!usedPositions[pos.y()][pos.x()]) {
				child[i] = pos;
				usedPositions[pos.y()][pos.x()] = true;
			}
		}
	}

	// For any unplaced genes, find available positions
	std::vector<Coordinate> availablePositions;
	for (size_t y = 0; y < field.rows; ++y) {
		for (size_t x = 0; x < field.cols; ++x) {
			if (!usedPositions[y][x]) {
				availablePositions.emplace_back(x, y);
			}
		}
	}

	// Shuffle available positions
	std::shuffle(availablePositions.begin(), availablePositions.end(), gen);

	// Place unplaced genes
	size_t availableIndex = 0;
	for (size_t i = 0; i < child.size(); ++i) {
		if (!child[i].isValid()) {
			if (availableIndex < availablePositions.size()) {
				child[i] = availablePositions[availableIndex++];
			}
			else {
				// No more available positions - this shouldn't happen with valid inputs
				// but we'll handle it gracefully
				child[i] = Coordinate::invalid();
			}
		}
	}

	// Add fillers if allowed
	if (field.allowFillers) {
		child.reserveFillers(availablePositions.size() - availableIndex);

		for (size_t i = availableIndex; i < availablePositions.size(); ++i) {
			child.addFiller(availablePositions[i]);
		}
	}

	return child;
}

// Implementation of the mutation function
inline void GeneticAlgorithm::mutate(Chromosome& chromosome, const Scheme& scheme) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0.0, 1.0);
	std::uniform_int_distribution<size_t> indexDist(0, chromosome.size() - 1);

	// Calculate gene mutation probability
	double geneMutationProb = 1.0 / chromosome.size();

	// For each gene, decide whether to mutate
	for (size_t i = 0; i < chromosome.size(); ++i) {
		if (dist(gen) < geneMutationProb) {
			// Choose another gene to swap with
			size_t j = indexDist(gen);

			// Don't swap a gene with itself
			if (i != j) {
				// Swap positions
				std::swap(chromosome[i], chromosome[j]);
			}
		}
	}

	// Optionally, also swap with fillers
	if (!chromosome.getFillers().empty() && dist(gen) < 0.2) {
		size_t geneIndex = indexDist(gen);
		std::uniform_int_distribution<size_t> fillerDist(0, chromosome.getFillers().size() - 1);
		size_t fillerIndex = fillerDist(gen);

		// Swap a gene with a filler
		std::swap(chromosome[geneIndex], chromosome.getFillers()[fillerIndex]);
	}
}

// Implementation of configureForType method
inline void SpecializedPopulation::configureForType(SpecializedPopType type) {
	MultiParametricFitness::Weights weights;

	// Reset all weights
	weights.wireLength = 0.0;
	weights.thermal = 0.0;
	weights.drc = 0.0;
	weights.power = 0.0;
	weights.parasitic = 0.0;

	// Set specific weights based on type
	switch (type) {
	case SpecializedPopType::WireLength:
		weights.wireLength = 1.0;
		break;
	case SpecializedPopType::Thermal:
		weights.thermal = 1.0;
		break;
	case SpecializedPopType::DRC:
		weights.drc = 1.0;
		break;
	case SpecializedPopType::Power:
		weights.power = 1.0;
		break;
	case SpecializedPopType::Parasitic:
		weights.parasitic = 1.0;
		break;
	case SpecializedPopType::WireThermal:
		weights.wireLength = 0.7;
		weights.thermal = 0.3;
		break;
	case SpecializedPopType::PowerDRC:
		weights.power = 0.4;
		weights.drc = 0.6;
		break;
	case SpecializedPopType::Full:
		weights.wireLength = 1.0;
		weights.thermal = 0.5;
		weights.drc = 2.0;
		weights.power = 0.3;
		weights.parasitic = 0.2;
		break;
	}

	m_fitness.setWeights(weights);
}

// Implementation of setFitnessParameters method
inline void SpecializedPopulation::setFitnessParameters(
	const MultiParametricFitness::ThermalModelParams& thermalParams,
	const MultiParametricFitness::PowerModelParams& powerParams,
	DRCChecker* drcChecker) {

	m_fitness.setThermalParams(thermalParams);
	m_fitness.setPowerParams(powerParams);
	m_fitness.setDRCChecker(drcChecker);
}

// Implementation of initialize method
inline void SpecializedPopulation::initialize(const Scheme& scheme) {
	// Fill with random placements
	m_population.fillWithRandomPlacements(scheme,
		[this](const Chromosome& chromosome, const Scheme& scheme) {
			return m_fitness.calculate(chromosome, scheme);
		});
}

// Implementation of evolve method
inline void SpecializedPopulation::evolve(
	const Scheme& scheme,
	const std::function<Chromosome(const Chromosome&, const Chromosome&, const Scheme&)>& crossoverFunc,
	const std::function<void(Chromosome&, const Scheme&)>& mutateFunc,
	double mutationProbability) {

	for (size_t i = 0; i < m_evolutions; i++) {
		// Filter out worst individuals
		size_t missingCount = m_population.filterByFitness(0.4); // Keep top 60%

		// Select parents for reproduction
		auto parentPairs = m_population.selectParentPairs();

		// Create offspring
		for (const auto& [parent1Idx, parent2Idx] : parentPairs) {
			const auto& parent1 = m_population[parent1Idx];
			const auto& parent2 = m_population[parent2Idx];

			// Create child through crossover
			Chromosome child = crossoverFunc(parent1.chromosome, parent2.chromosome, scheme);

			// Potentially mutate
			if (std::rand() % 100 < mutationProbability * 100) {
				mutateFunc(child, scheme);
			}

			// Calculate fitness
			double fitness = m_fitness.calculate(child, scheme);

			// Add to population
			m_population.addIndividual(std::move(child), fitness);
		}
	}
}

// Implementation of getBestIndividuals method
inline std::vector<Individual> SpecializedPopulation::getBestIndividuals(double percentage) const {
	// Calculate how many to take
	size_t count = static_cast<size_t>(m_population.size() * percentage);
	if (count == 0) count = 1; // Take at least one

	// Create a vector for indices
	std::vector<size_t> indices(m_population.size());
	for (size_t i = 0; i < indices.size(); ++i) {
		indices[i] = i;
	}

	// Sort indices by fitness
	std::sort(indices.begin(), indices.end(),
		[this](size_t a, size_t b) {
			return m_population[a].fitness < m_population[b].fitness;
		});

	// Copy best individuals
	std::vector<Individual> best;
	best.reserve(count);

	for (size_t i = 0; i < count && i < indices.size(); ++i) {
		size_t idx = indices[i];
		// Create a deep copy of the chromosome
		Chromosome chromosomeCopy(m_population[idx].chromosome);
		best.emplace_back(std::move(chromosomeCopy), m_population[idx].fitness);
	}

	return best;
}

// Implementation of getRandomIndividuals method
inline std::vector<Individual> SpecializedPopulation::getRandomIndividuals(double percentage) const {
	// Calculate how many to take
	size_t count = static_cast<size_t>(m_population.size() * percentage);
	if (count == 0) count = 1; // Take at least one

	// Create a vector of indices
	std::vector<size_t> indices(m_population.size());
	for (size_t i = 0; i < indices.size(); i++) {
		indices[i] = i;
	}

	// Shuffle the indices
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);

	// Copy random individuals
	std::vector<Individual> random;
	random.reserve(count);

	for (size_t i = 0; i < count && i < indices.size(); i++) {
		size_t idx = indices[i];
		// Create a deep copy of the chromosome
		Chromosome chromosomeCopy(m_population[idx].chromosome);
		random.emplace_back(std::move(chromosomeCopy), m_population[idx].fitness);
	}

	return random;
}

// Implementation of setFitnessParameters method for manager
inline void SpecializedPopulationsManager::setFitnessParameters(
	const MultiParametricFitness::ThermalModelParams& thermalParams,
	const MultiParametricFitness::PowerModelParams& powerParams,
	DRCChecker* drcChecker) {

	for (auto& pop : m_populations) {
		pop.setFitnessParameters(thermalParams, powerParams, drcChecker);
	}
}

// Implementation of initializeAll method
inline void SpecializedPopulationsManager::initializeAll(const Scheme& scheme) {
	for (auto& pop : m_populations) {
		pop.initialize(scheme);
	}
}

// Implementation of evolveAll method
inline void SpecializedPopulationsManager::evolveAll(
	const Scheme& scheme,
	const std::function<Chromosome(const Chromosome&, const Chromosome&, const Scheme&)>& crossoverFunc,
	const std::function<void(Chromosome&, const Scheme&)>& mutateFunc,
	double mutationProbability) {

	std::vector<std::thread> threads;
	threads.reserve(m_populations.size());

	// Launch a thread for each population
	for (auto& pop : m_populations) {
		threads.emplace_back([&pop, &scheme, &crossoverFunc, &mutateFunc, mutationProbability]() {
			pop.evolve(scheme, crossoverFunc, mutateFunc, mutationProbability);
			});
	}

	// Wait for all threads to finish
	for (auto& thread : threads) {
		thread.join();
	}
}

// Implementation of createMergedPopulation method
inline Population SpecializedPopulationsManager::createMergedPopulation(
	const Scheme& scheme,
	const MultiParametricFitness& fitness,
	double bestPercentage,
	double randomPercentage) {

	// Create a new population
	Population mergedPop(m_basePopulationSize);

	// Get best individuals from each population
	for (auto& pop : m_populations) {
		auto best = pop.getBestIndividuals(bestPercentage);
		auto random = pop.getRandomIndividuals(randomPercentage);

		// Add best individuals
		for (auto& ind : best) {
			if (mergedPop.size() < mergedPop.capacity()) {
				// Recalculate fitness with the full fitness function
				double newFitness = fitness.calculate(ind.chromosome, scheme);

				// Create a new chromosome (deep copy) and move it to the population
				Chromosome newChromosome(ind.chromosome);
				mergedPop.addIndividual(std::move(newChromosome), newFitness);
			}
		}

		// Add random individuals
		for (auto& ind : random) {
			if (mergedPop.size() < mergedPop.capacity()) {
				// Recalculate fitness with the full fitness function
				double newFitness = fitness.calculate(ind.chromosome, scheme);

				// Create a new chromosome (deep copy) and move it to the population
				Chromosome newChromosome(ind.chromosome);
				mergedPop.addIndividual(std::move(newChromosome), newFitness);
			}
		}
	}

	return mergedPop;
}

// Implementation of runWithSpecializedPopulations method
inline Chromosome GeneticAlgorithm::runWithSpecializedPopulations(const Scheme& scheme) {
	auto startTime = std::chrono::high_resolution_clock::now();

	// Create specialized populations manager
	SpecializedPopulationsManager manager(m_params.populationSize);

	// Add specialized populations
	manager.addSpecializedPopulation(SpecializedPopType::WireLength);
	manager.addSpecializedPopulation(SpecializedPopType::Thermal);
	manager.addSpecializedPopulation(SpecializedPopType::DRC);
	manager.addSpecializedPopulation(SpecializedPopType::WireThermal);
	manager.addSpecializedPopulation(SpecializedPopType::PowerDRC);

	// Configure specialized populations
	manager.setEvolutions(m_params.specializedEvolutions);
	manager.setFitnessParameters(
		m_fitness.getThermalParams(),
		m_fitness.getPowerParams(),
		&m_drcChecker
	);

	// Initialize and evolve specialized populations
	if (m_verbose) {
		std::cout << "Initializing specialized populations...\n";
	}
	manager.initializeAll(scheme);

	if (m_verbose) {
		std::cout << "Evolving specialized populations...\n";
	}
	manager.evolveAll(scheme,
		[](const Chromosome& p1, const Chromosome& p2, const Scheme& s) {
			return crossover(p1, p2, s);
		},
		[](Chromosome& c, const Scheme& s) {
			mutate(c, s);
		},
		m_params.mutationProbability);

	// Create merged population
	if (m_verbose) {
		std::cout << "Merging specialized populations...\n";
	}
	Population population = manager.createMergedPopulation(
		scheme,
		m_fitness,
		m_params.specializedBestPercentage,
		m_params.specializedRandomPercentage
	);

	// Fill remaining spots with random placements if needed
	population.fillWithRandomPlacements(scheme,
		[this](const Chromosome& chromosome, const Scheme& scheme) {
			return m_fitness.calculate(chromosome, scheme);
		});

	// Record initial fitness
	m_stats.initialFitness = population.calculateAverageFitness();

	// Continue with standard evolution
	return runStandardEvolution(population, scheme);
}

// Implementation of runStandardEvolution method
inline Chromosome GeneticAlgorithm::runStandardEvolution(Population& population, const Scheme& scheme) {
	auto startTime = std::chrono::high_resolution_clock::now();

	// If population is empty, initialize it
	if (population.size() == 0) {
		if (m_verbose) {
			std::cout << "Initializing population...\n";
		}
		population = Population(m_params.populationSize);
		population.fillWithRandomPlacements(scheme,
			[this](const Chromosome& chromosome, const Scheme& scheme) {
				return m_fitness.calculate(chromosome, scheme);
			});

		// Record initial fitness
		m_stats.initialFitness = population.calculateAverageFitness();
	}

	double prevAvgFitness = m_stats.initialFitness;
	double iterationsWithoutImprovement = 0;

	// Main evolution loop
	for (size_t i = 0; i < m_params.maxIterations; ++i) {
		// Filter population
		size_t missingCount = population.filterByFitness(m_params.filterCoefficient);

		// Check if we have enough individuals
		if (population.size() < 2) {
			// Refill with random placements to avoid degeneration
			population.fillWithRandomPlacements(scheme,
				[this](const Chromosome& chromosome, const Scheme& scheme) {
					return m_fitness.calculate(chromosome, scheme);
				});
		}

		// Select parents
		auto parentPairs = population.selectParentPairs();

		// Create offspring
		for (const auto& [parent1Idx, parent2Idx] : parentPairs) {
			// Безопасная проверка индексов
			if (parent1Idx >= population.size() || parent2Idx >= population.size()) {
				continue;
			}

			const auto& parent1 = population[parent1Idx];
			const auto& parent2 = population[parent2Idx];

			// Create child through crossover
			Chromosome child = crossover(parent1.chromosome, parent2.chromosome, scheme);

			// Potentially mutate
			if (std::rand() % 100 < m_params.mutationProbability * 100) {
				mutate(child, scheme);
			}

			// Calculate fitness
			double fitness = m_fitness.calculate(child, scheme);

			// Add to population
			population.addIndividual(std::move(child), fitness);
		}

		// Periodically check DRC and filter worst violations
		if (i % m_params.drcCheckInterval == 0 && i > 0) {
			population = m_drcChecker.filterPopulationByDRC(
				population,
				scheme,
				m_params.drcFilterPercentage
			);

			// Fill back up to capacity
			population.fillWithRandomPlacements(scheme,
				[this](const Chromosome& chromosome, const Scheme& scheme) {
					return m_fitness.calculate(chromosome, scheme);
				});
		}

		// Calculate current average fitness
		double currentAvgFitness = population.calculateAverageFitness();
		//std::cout << "Iteration: " << i << ", Average Fitness: " << currentAvgFitness << std::endl;

		// Check for convergence
		double relativeImprovement = 0.0;
		if (prevAvgFitness > 0) {
			relativeImprovement = std::abs(prevAvgFitness - currentAvgFitness) / prevAvgFitness;
		}

		if (relativeImprovement < m_params.convergenceThreshold) {
			iterationsWithoutImprovement++;

			// Stop if no improvement for a while
			if (iterationsWithoutImprovement >= 10) {
				m_stats.iterations = i + 1;
				break;
			}
		}
		else {
			iterationsWithoutImprovement = 0;
		}

		prevAvgFitness = currentAvgFitness;

		// Print progress
		if (m_verbose && (i % 10 == 0 || i == m_params.maxIterations - 1)) {
			population.sortByFitness();
			if (!population.size()) {
				std::cout << "Warning: Population is empty on iteration " << i << std::endl;
				continue;
			}
			std::cout << "Iteration " << i <<
				", Avg Fitness = " << currentAvgFitness <<
				", Best = " << population[0].fitness << "\n";
		}

		// Safety check - prevent infinite loops
		if (i > 0 && i % 50 == 0) {
			// Inject random individuals to maintain diversity
			size_t randomCount = static_cast<size_t>(population.capacity() * 0.1); // 10% new random individuals
			for (size_t j = 0; j < randomCount && population.size() < population.capacity(); j++) {
				Chromosome randomChromosome(scheme);
				randomChromosome.generateRandomPlacement(scheme);
				double fitness = m_fitness.calculate(randomChromosome, scheme);
				population.addIndividual(std::move(randomChromosome), fitness);
			}
		}
	}

	// Final sort to get the best individual
	population.sortByFitness();

	// Safety check - make sure we have at least one individual
	if (population.size() == 0) {
		if (m_verbose) {
			std::cout << "Warning: No individuals in final population. Generating a random solution.\n";
		}

		// Generate a random solution
		Chromosome randomSolution(scheme);
		randomSolution.generateRandomPlacement(scheme);
		population.addIndividual(std::move(randomSolution),
			m_fitness.calculate(population[0].chromosome, scheme));
	}

	// Record final statistics
	m_stats.finalFitness = population.calculateAverageFitness();

	// Calculate individual fitness components for the best chromosome
	const Chromosome& bestChromosome = population[0].chromosome;
	m_stats.bestWireLengthFitness = m_fitness.calculateWireLengthFitness(bestChromosome, scheme);
	m_stats.bestThermalFitness = m_fitness.calculateThermalFitness(bestChromosome, scheme);
	m_stats.bestDRCFitness = m_fitness.calculateDRCFitness(bestChromosome, scheme);
	m_stats.bestPowerFitness = m_fitness.calculatePowerFitness(bestChromosome, scheme);
	m_stats.bestParasiticFitness = m_fitness.calculateParasiticFitness(bestChromosome, scheme);

	// Record execution time
	auto endTime = std::chrono::high_resolution_clock::now();
	m_stats.executionTime = std::chrono::duration_cast<std::chrono::milliseconds>(
		endTime - startTime
	);

	return bestChromosome;
}