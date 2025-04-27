#pragma once

#include "Chromosome.h"
#include "Scheme.h"
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <cassert>

// Structure to represent an individual in the population
struct Individual {
	Chromosome chromosome;
	double fitness;

	// Default constructor
	Individual() : fitness(0.0) {}

	// Constructor with chromosome and fitness
	Individual(Chromosome&& chrom, double fit)
		: chromosome(std::move(chrom)), fitness(fit) {
	}

	// Copy constructor - создает глубокую копию хромосомы
	Individual(const Individual& other)
		: chromosome(other.chromosome), fitness(other.fitness) {
	}

	// Move constructor
	Individual(Individual&& other) = default;

	// Copy assignment
	Individual& operator=(const Individual& other) {
		if (this != &other) {
			chromosome = other.chromosome;
			fitness = other.fitness;
		}
		return *this;
	}

	// Move assignment
	Individual& operator=(Individual&& other) = default;
};

// Alias for a collection of individuals
using PopulationVector = std::vector<Individual>;

// Class representing a population of chromosomes
class Population {
private:
	PopulationVector m_individuals;
	size_t m_capacity;

public:
	// Constructor
	Population(size_t capacity)
		: m_capacity(capacity) {
		m_individuals.reserve(capacity);
	}

	// Access operators
	Individual& operator[](size_t index) {
		return m_individuals[index];
	}

	const Individual& operator[](size_t index) const {
		return m_individuals[index];
	}

	// Get current size
	size_t size() const {
		return m_individuals.size();
	}

	// Get capacity
	size_t capacity() const {
		return m_capacity;
	}

	// Get number of individuals missing from capacity
	size_t missingCount() const {
		return capacity() - size();
	}

	// Add an individual to the population
	void addIndividual(Chromosome&& chromosome, double fitness) {
		m_individuals.emplace_back(std::move(chromosome), fitness);
	}

	// Fill population with random placements
	void fillWithRandomPlacements(const Scheme& scheme,
		const std::function<double(const Chromosome&, const Scheme&)>& fitnessFunction) {
		while (size() < capacity()) {
			Chromosome chromosome(scheme);
			chromosome.generateRandomPlacement(scheme);

			// Calculate fitness
			double fitness = fitnessFunction(chromosome, scheme);

			// Add to population
			addIndividual(std::move(chromosome), fitness);
		}
	}

	// Filter the population by fitness (keep top individuals)
	size_t filterByFitness(double filtrationCoeff) {
		assert(filtrationCoeff >= 0.0 && filtrationCoeff <= 1.0);

		// Calculate new population size
		size_t newSize = static_cast<size_t>((1.0 - filtrationCoeff) * size());
		if (newSize == 0) newSize = 1; // Keep at least one individual

		// Sort population by fitness (lower is better)
		sortByFitness();

		// Remove individuals with higher fitness values
		m_individuals.resize(newSize);

		return missingCount();
	}

	// Sort population by fitness (lower is better)
	void sortByFitness() {
		std::sort(m_individuals.begin(), m_individuals.end(),
			[](const Individual& a, const Individual& b) {
				return a.fitness < b.fitness;
			});
	}

	// Calculate probabilities for parent selection
	std::vector<double> getParentSelectionProbabilities() const {
		// Calculate sum of fitness values
		double fitnessSum = 0.0;
		for (const auto& individual : m_individuals) {
			fitnessSum += individual.fitness;
		}

		// Calculate probabilities (lower fitness = higher probability)
		std::vector<double> probabilities(size(), 0.0);

		// Check if fitnessSum is valid and there are valid individuals
		if (fitnessSum <= 0.0 || size() <= 1) {
			// Use uniform distribution if fitness sum is zero or there's only one individual
			double uniformProb = 1.0 / size();
			std::fill(probabilities.begin(), probabilities.end(), uniformProb);
			return probabilities;
		}

		for (size_t i = 0; i < size(); ++i) {
			// Invert probability so that lower fitness results in higher probability
			probabilities[i] = (fitnessSum - m_individuals[i].fitness) / (fitnessSum * (size() - 1));

			// Make sure probability is positive
			probabilities[i] = std::max(0.000001, probabilities[i]);
		}

		return probabilities;
	}

	// Select parent pairs for crossover
	std::vector<std::pair<size_t, size_t>> selectParentPairs() const {
		auto probabilities = getParentSelectionProbabilities();

		// Make sure we have valid probabilities
		double probSum = 0.0;
		for (auto p : probabilities) {
			probSum += p;
		}

		if (probSum <= 0.0) {
			// If still no valid probabilities, use uniform
			std::fill(probabilities.begin(), probabilities.end(), 1.0 / probabilities.size());
		}

		std::vector<std::pair<size_t, size_t>> parentPairs;
		parentPairs.reserve(missingCount());

		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<> dist(probabilities.begin(), probabilities.end());

		for (size_t i = 0; i < missingCount(); ++i) {
			size_t parent1 = dist(gen);
			size_t parent2;

			// If we only have one individual, just use it twice
			if (size() <= 1) {
				parent2 = parent1;
			}
			else {
				// Ensure parent2 is different from parent1
				do {
					parent2 = dist(gen);
				} while (parent2 == parent1);
			}

			parentPairs.emplace_back(parent1, parent2);
		}

		return parentPairs;
	}

	// Calculate total fitness
	double calculateTotalFitness() const {
		double total = 0.0;
		for (const auto& individual : m_individuals) {
			total += individual.fitness;
		}
		return total;
	}

	// Calculate average fitness
	double calculateAverageFitness() const {
		return size() > 0 ? calculateTotalFitness() / size() : 0.0;
	}

	// Get best (lowest fitness) individual
	const Individual& getBestIndividual() const {
		assert(!m_individuals.empty());

		auto it = std::min_element(m_individuals.begin(), m_individuals.end(),
			[](const Individual& a, const Individual& b) {
				return a.fitness < b.fitness;
			});

		return *it;
	}

	// Get worst (highest fitness) individual
	const Individual& getWorstIndividual() const {
		assert(!m_individuals.empty());

		auto it = std::max_element(m_individuals.begin(), m_individuals.end(),
			[](const Individual& a, const Individual& b) {
				return a.fitness < b.fitness;
			});

		return *it;
	}

	// Clear the population
	void clear() {
		m_individuals.clear();
	}
};