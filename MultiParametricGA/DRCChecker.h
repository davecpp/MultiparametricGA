#pragma once

#include "Cell.h"
#include "Chromosome.h"
#include "Coordinate.h"
#include "Scheme.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <iostream>

// Class for Design Rule Checking (DRC)
class DRCChecker {
public:
	// DRC parameters
	struct DRCParams {
		int32_t minDistance = 2;
		double overlapPenalty = 10.0;
		double distancePenalty = 5.0;
	};

private:
	DRCParams m_params;

public:
	// Constructor
	DRCChecker() = default;

	// Constructor with parameters
	DRCChecker(const DRCParams& params) : m_params(params) {}

	// Set parameters
	void setParams(const DRCParams& params) {
		m_params = params;
	}

	// Get parameters
	const DRCParams& getParams() const {
		return m_params;
	}

	// Calculate DRC penalty for a chromosome
	double calculatePenalty(const Chromosome& chromosome, const Scheme& scheme) const {
		try {
			const auto& cells = scheme.getCells();
			double totalPenalty = 0.0;

			// Check for overlaps and minimum distance violations
			for (size_t i = 0; i < cells.size() && i < chromosome.size(); i++) {
				const Coordinate& pos1 = chromosome[i];
				if (!pos1.isValid()) continue;

				const Cell& cell1 = cells[i];
				// Получаем размеры ячейки через getSize()
				auto [width1, height1] = cell1.getSize();

				for (size_t j = i + 1; j < cells.size() && j < chromosome.size(); j++) {
					const Coordinate& pos2 = chromosome[j];
					if (!pos2.isValid()) continue;

					const Cell& cell2 = cells[j];
					// Получаем размеры ячейки через getSize()
					auto [width2, height2] = cell2.getSize();

					// Calculate distance between cells
					int32_t distance = Coordinate::distance(pos1, pos2);

					// Check for overlaps
					bool overlapX = (pos1.x() < pos2.x() + width2) && (pos2.x() < pos1.x() + width1);
					bool overlapY = (pos1.y() < pos2.y() + height2) && (pos2.y() < pos1.y() + height1);

					if (overlapX && overlapY) {
						// Calculate overlap area
						int32_t overlapWidth = std::min(pos1.x() + width1, pos2.x() + width2) -
							std::max(pos1.x(), pos2.x());
						int32_t overlapHeight = std::min(pos1.y() + height1, pos2.y() + height2) -
							std::max(pos1.y(), pos2.y());

						double overlapArea = overlapWidth * overlapHeight;
						totalPenalty += m_params.overlapPenalty * overlapArea;
					}
					// Check minimum distance violation
					else if (distance < m_params.minDistance) {
						totalPenalty += m_params.distancePenalty * (m_params.minDistance - distance);
					}
				}
			}

			return totalPenalty;
		}
		catch (const std::exception& e) {
			std::cerr << "Error in DRCChecker::calculatePenalty: " << e.what() << std::endl;
			return 0.0;
		}
	}

	// Filter population by DRC violations
	Population filterPopulationByDRC(const Population& population, const Scheme& scheme, double filterPercentage) const {
		try {
			// Calculate DRC penalties for all individuals
			std::vector<std::pair<size_t, double>> penalties;
			for (size_t i = 0; i < population.size(); i++) {
				double penalty = calculatePenalty(population[i].chromosome, scheme);
				penalties.push_back({ i, penalty });
			}

			// Sort by penalty (highest first)
			std::sort(penalties.begin(), penalties.end(),
				[](const auto& a, const auto& b) {
					return a.second > b.second;
				});

			// Determine how many to filter out
			size_t filterCount = static_cast<size_t>(population.size() * filterPercentage);
			if (filterCount == 0 && !penalties.empty() && penalties[0].second > 0) {
				filterCount = 1; // Filter at least one if there are DRC violations
			}

			// Create a set of indices to filter out
			std::unordered_set<size_t> indicesToFilter;
			for (size_t i = 0; i < filterCount && i < penalties.size(); i++) {
				indicesToFilter.insert(penalties[i].first);
			}

			// Create new population without filtered individuals
			Population filteredPopulation(population.capacity());
			for (size_t i = 0; i < population.size(); i++) {
				if (indicesToFilter.find(i) == indicesToFilter.end()) {
					// Create a deep copy of the chromosome
					Chromosome chromosomeCopy(population[i].chromosome);
					filteredPopulation.addIndividual(
						std::move(chromosomeCopy),
						population[i].fitness
					);
				}
			}

			return filteredPopulation;
		}
		catch (const std::exception& e) {
			std::cerr << "Error in DRCChecker::filterPopulationByDRC: " << e.what() << std::endl;
			return population; // Return the original population in case of error
		}
	}
};