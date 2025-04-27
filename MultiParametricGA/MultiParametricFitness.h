#pragma once

#include "Cell.h"
#include "Chromosome.h"
#include "Coordinate.h"
#include "Scheme.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

// Forward declaration
class DRCChecker;

// Class implementing multi-parametric fitness function
class MultiParametricFitness {
public:
	// Weight structure for different components
	struct Weights {
		double wireLength = 1.0;
		double thermal = 0.5;
		double drc = 2.0;
		double power = 0.3;
		double parasitic = 0.2;
	};

	// Thermal model parameters
	struct ThermalModelParams {
		double maxAllowedTemp = 100.0;
		double ambientTemp = 25.0;
		double thermalConductivity = 150.0;
		double powerDensity = 0.1;
	};

	// Power model parameters
	struct PowerModelParams {
		double cellPower = 0.1;
		double wireCapacitance = 0.02;
		double voltageFactor = 1.0;
		double frequency = 1.0;
		double leakageFactor = 0.05;
	};

private:
	Weights m_weights;
	ThermalModelParams m_thermalParams;
	PowerModelParams m_powerParams;
	DRCChecker* m_drcChecker = nullptr;

public:
	// Constructor
	MultiParametricFitness() = default;

	// Constructor with parameters
	MultiParametricFitness(const Weights& weights,
		const ThermalModelParams& thermalParams,
		const PowerModelParams& powerParams)
		: m_weights(weights),
		m_thermalParams(thermalParams),
		m_powerParams(powerParams) {
	}

	// Set weights
	void setWeights(const Weights& weights) {
		m_weights = weights;
	}

	// Get weights
	const Weights& getWeights() const {
		return m_weights;
	}

	// Set thermal parameters
	void setThermalParams(const ThermalModelParams& params) {
		m_thermalParams = params;
	}

	// Get thermal parameters
	const ThermalModelParams& getThermalParams() const {
		return m_thermalParams;
	}

	// Set power parameters
	void setPowerParams(const PowerModelParams& params) {
		m_powerParams = params;
	}

	// Get power parameters
	const PowerModelParams& getPowerParams() const {
		return m_powerParams;
	}

	// Set DRC checker
	void setDRCChecker(DRCChecker* checker) {
		m_drcChecker = checker;
	}

	// Get DRC checker
	DRCChecker* getDRCChecker() const {
		return m_drcChecker;
	}

	// Calculate wire length component
	double calculateWireLengthFitness(const Chromosome& chromosome, const Scheme& scheme) const {
		try {
			const auto& connections = scheme.getConnections();

			// Safety checks
			if (connections.rows() == 0 || connections.cols() == 0 || chromosome.size() == 0) {
				return 0.0;
			}

			double wireLength = 0.0;
			const size_t cellsCount = std::min(chromosome.size(), connections.rows());

			for (size_t i = 0; i < cellsCount; i++) {
				for (size_t j = 0; j < cellsCount; j++) {
					double weight = connections[i][j];
					if (weight > 0) {
						int32_t distance = Coordinate::distance(chromosome[i], chromosome[j]);
						wireLength += weight * distance;
					}
				}
			}

			return wireLength / cellsCount;
		}
		catch (const std::exception& e) {
			std::cerr << "Error in calculateWireLengthFitness: " << e.what() << std::endl;
			return 0.0;
		}
	}

	// Calculate thermal component
	double calculateThermalFitness(const Chromosome& chromosome, const Scheme& scheme) const {
		try {
			const auto& cells = scheme.getCells();
			const auto& field = scheme.getField();

			// Safety check
			if (field.rows == 0 || field.cols == 0 || chromosome.size() == 0) {
				return 0.0;
			}

			// Create a thermal grid
			Matrix<double> thermalGrid(field.rows, field.cols, m_thermalParams.ambientTemp);

			// Calculate thermal distribution
			for (size_t i = 0; i < std::min(cells.size(), chromosome.size()); i++) {
				const Coordinate& pos = chromosome[i];
				if (!pos.isValid()) continue;

				// Safety check for coordinates
				if (pos.x() < 0 || pos.y() < 0 ||
					pos.x() >= static_cast<int>(field.cols) ||
					pos.y() >= static_cast<int>(field.rows)) {
					continue;
				}

				const Cell& cell = cells[i];
				double cellPower = std::max(0.001, cell.getThermalValue()) * m_thermalParams.powerDensity;

				// Apply heat to surrounding area (simple diffusion model)
				for (int dy = -2; dy <= 2; dy++) {
					for (int dx = -2; dx <= 2; dx++) {
						int x = pos.x() + dx;
						int y = pos.y() + dy;

						if (x >= 0 && x < static_cast<int>(field.cols) &&
							y >= 0 && y < static_cast<int>(field.rows)) {

							// Heat diminishes with distance (inverse square)
							double distance = std::sqrt(dx * dx + dy * dy);
							double heatContribution = cellPower / (1.0 + distance * distance);

							thermalGrid[y][x] += heatContribution;
						}
					}
				}
			}

			// Calculate thermal metrics
			double maxTemp = m_thermalParams.ambientTemp;
			double sumTemp = 0.0;
			double sumSquaredTemp = 0.0;

			for (size_t y = 0; y < field.rows; y++) {
				for (size_t x = 0; x < field.cols; x++) {
					double temp = thermalGrid[y][x];
					maxTemp = std::max(maxTemp, temp);
					sumTemp += temp;
					sumSquaredTemp += temp * temp;
				}
			}

			double avgTemp = sumTemp / (field.rows * field.cols);
			double tempVariance = (sumSquaredTemp / (field.rows * field.cols)) - (avgTemp * avgTemp);

			// Higher max temperature and variance = worse placement
			return (maxTemp - m_thermalParams.ambientTemp) + std::sqrt(std::max(0.0, tempVariance));
		}
		catch (const std::exception& e) {
			std::cerr << "Error in calculateThermalFitness: " << e.what() << std::endl;
			return 0.0;
		}
	}

	// Calculate power consumption component
	double calculatePowerFitness(const Chromosome& chromosome, const Scheme& scheme) const {
		try {
			const auto& cells = scheme.getCells();
			const auto& connections = scheme.getConnections();

			// Safety check
			if (connections.rows() == 0 || connections.cols() == 0 ||
				chromosome.size() == 0 || cells.size() == 0) {
				return 0.0;
			}

			double dynamicPower = 0.0;
			double leakagePower = 0.0;

			// Calculate dynamic power based on connections and wire lengths
			for (size_t i = 0; i < std::min(cells.size(), chromosome.size()); i++) {
				if (!chromosome[i].isValid()) continue;

				// Base cell power (защита от отрицательных значений)
				dynamicPower += m_powerParams.cellPower * std::max(0.001, cells[i].getPowerDensity());

				// Power from interconnections
				for (size_t j = 0; j < std::min(cells.size(), chromosome.size()); j++) {
					if (!chromosome[j].isValid()) continue;

					if (i < connections.rows() && j < connections.cols()) {
						double connectionWeight = connections[i][j];
						if (connectionWeight > 0) {
							int32_t distance = Coordinate::distance(chromosome[i], chromosome[j]);
							double wireCap = distance * m_powerParams.wireCapacitance;

							// P = α·C·V²·f
							double activityFactor = connectionWeight;
							dynamicPower += activityFactor * wireCap *
								std::pow(m_powerParams.voltageFactor, 2) *
								m_powerParams.frequency;
						}
					}
				}

				// Simplified leakage model based on thermal value
				double cellTemp = std::max(0.001, cells[i].getThermalValue()) + m_thermalParams.ambientTemp;
				leakagePower += m_powerParams.leakageFactor * cells[i].getPowerDensity() *
					std::exp((cellTemp - m_thermalParams.ambientTemp) / 10.0);
			}

			return dynamicPower + leakagePower;
		}
		catch (const std::exception& e) {
			std::cerr << "Error in calculatePowerFitness: " << e.what() << std::endl;
			return 0.0;
		}
	}

	// Calculate parasitic effects component
	double calculateParasiticFitness(const Chromosome& chromosome, const Scheme& scheme) const {
		try {
			const auto& cells = scheme.getCells();
			const auto& connections = scheme.getConnections();

			// Safety check
			if (connections.rows() == 0 || connections.cols() == 0 ||
				chromosome.size() == 0 || cells.size() == 0) {
				return 0.0;
			}

			double parasiticEffects = 0.0;

			// Parasitic effects increase with wire length and decrease with distance between wires
			for (size_t i = 0; i < std::min(cells.size(), chromosome.size()); i++) {
				if (!chromosome[i].isValid()) continue;

				for (size_t j = i + 1; j < std::min(cells.size(), chromosome.size()); j++) {
					if (!chromosome[j].isValid()) continue;

					if (i < connections.rows() && j < connections.cols() && connections[i][j] > 0) {
						int32_t wireLength = Coordinate::distance(chromosome[i], chromosome[j]);

						// Check for parallel wires (greater coupling)
						for (size_t k = 0; k < std::min(cells.size(), chromosome.size()); k++) {
							if (k == i || k == j || !chromosome[k].isValid()) continue;

							if (i < connections.rows() && k < connections.cols() && connections[i][k] > 0) {
								int32_t parallelDistance = std::abs(
									Coordinate::distance(chromosome[i], chromosome[k]) -
									Coordinate::distance(chromosome[j], chromosome[k])
								);

								// Closer parallel wires = more coupling
								if (parallelDistance < 3) {
									parasiticEffects += wireLength * (3 - parallelDistance) / 3.0;
								}
							}
						}
					}
				}
			}

			return parasiticEffects;
		}
		catch (const std::exception& e) {
			std::cerr << "Error in calculateParasiticFitness: " << e.what() << std::endl;
			return 0.0;
		}
	}

	// Calculate DRC fitness component
	double calculateDRCFitness(const Chromosome& chromosome, const Scheme& scheme) const {
		try {
			if (m_drcChecker) {
				return m_drcChecker->calculatePenalty(chromosome, scheme);
			}
			return 0.0;
		}
		catch (const std::exception& e) {
			std::cerr << "Error in calculateDRCFitness: " << e.what() << std::endl;
			return 0.0;
		}
	}

	// Calculate complete fitness value (lower is better)
	double calculate(const Chromosome& chromosome, const Scheme& scheme) const {
		try {
			// Calculate individual components
			double wireLengthComponent = calculateWireLengthFitness(chromosome, scheme);
			double thermalComponent = calculateThermalFitness(chromosome, scheme);
			double powerComponent = calculatePowerFitness(chromosome, scheme);
			double parasiticComponent = calculateParasiticFitness(chromosome, scheme);

			// DRC component will be calculated separately if DRC checker is provided
			double drcComponent = 0.0;
			if (m_drcChecker) {
				drcComponent = calculateDRCFitness(chromosome, scheme);
			}

			// Combine with weights
			double fitness =
				m_weights.wireLength * wireLengthComponent +
				m_weights.thermal * thermalComponent +
				m_weights.drc * drcComponent +
				m_weights.power * powerComponent +
				m_weights.parasitic * parasiticComponent;

			return fitness;
		}
		catch (const std::exception& e) {
			std::cerr << "Error in calculate: " << e.what() << std::endl;
			return std::numeric_limits<double>::max(); // Возвращаем максимально плохое значение фитнеса
		}
	}
};