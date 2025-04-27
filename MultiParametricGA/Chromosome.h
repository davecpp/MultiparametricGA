#pragma once

#include "Cell.h"
#include "Coordinate.h"
#include "Scheme.h"
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <unordered_set>
#include <nlohmann/json.hpp>

// Class representing a chromosome (placement solution)
class Chromosome {
public:
    // Cell-Coordinate mapping types
    using CellPositions = std::vector<Coordinate>;
    using FillerPositions = std::vector<Coordinate>;

private:
    CellPositions m_positions;  // Cell positions (index = cell ID, value = position)
    FillerPositions m_fillers;  // Filler cell positions

public:
    // Default constructor
    Chromosome() = default;

    // Constructor with scheme initialization
    Chromosome(const Scheme& scheme)
        : m_positions(scheme.getCells().size(), Coordinate::invalid()) {
    }

    // Copy constructor
    Chromosome(const Chromosome& other)
        : m_positions(other.m_positions), m_fillers(other.m_fillers) {
    }

    // Move constructor
    Chromosome(Chromosome&& other) = default;

    // Assignment operators
    Chromosome& operator=(const Chromosome& other) = default;
    Chromosome& operator=(Chromosome&& other) = default;

    // Access operators
    Coordinate& operator[](size_t cellIndex) {
        return m_positions[cellIndex];
    }

    const Coordinate& operator[](size_t cellIndex) const {
        return m_positions[cellIndex];
    }

    // Get size (number of cells)
    size_t size() const {
        return m_positions.size();
    }

    // Get fillers
    const FillerPositions& getFillers() const {
        return m_fillers;
    }

    FillerPositions& getFillers() {
        return m_fillers;
    }

    // Reserve space for fillers
    void reserveFillers(size_t size) {
        m_fillers.reserve(size);
    }

    // Add a filler cell
    void addFiller(const Coordinate& position) {
        m_fillers.push_back(position);
    }

    // Generate a random placement
    void generateRandomPlacement(const Scheme& scheme) {
        const auto& field = scheme.getField();
        const auto& cells = scheme.getCells();

        // Create a list of all possible positions
        std::vector<Coordinate> availablePositions;
        for (size_t row = 0; row < field.rows; ++row) {
            for (size_t col = 0; col < field.cols; ++col) {
                availablePositions.emplace_back(col, row);
            }
        }

        // Shuffle the positions
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(availablePositions.begin(), availablePositions.end(), g);

        // Assign positions to cells
        for (size_t i = 0; i < cells.size(); ++i) {
            if (i < availablePositions.size()) {
                m_positions[i] = availablePositions[i];
            }
            else {
                // Not enough positions
                m_positions[i] = Coordinate::invalid();
            }
        }

        // If fillers are allowed, assign remaining positions to fillers
        if (field.allowFillers) {
            m_fillers.clear();
            for (size_t i = cells.size(); i < availablePositions.size(); ++i) {
                m_fillers.push_back(availablePositions[i]);
            }
        }
    }

    // Check if the chromosome is valid
    bool isValid(const Scheme& scheme) const {
        const auto& field = scheme.getField();

        // Check if all cells have valid positions
        for (const auto& pos : m_positions) {
            if (!pos.isValid()) {
                return false;
            }

            // Check if position is within field boundaries
            if (pos.x() >= static_cast<int32_t>(field.cols) ||
                pos.y() >= static_cast<int32_t>(field.rows) ||
                pos.x() < 0 || pos.y() < 0) {
                return false;
            }
        }

        // Check for duplicate positions
        std::unordered_set<uint64_t> positionSet;
        for (const auto& pos : m_positions) {
            // Create a unique hash for the position
            uint64_t hash = (static_cast<uint64_t>(pos.x()) << 32) | static_cast<uint64_t>(pos.y());

            // Check if this position is already used
            if (positionSet.find(hash) != positionSet.end()) {
                return false;
            }

            positionSet.insert(hash);
        }

        // Check filler positions
        for (const auto& pos : m_fillers) {
            // Check if position is within field boundaries
            if (pos.x() >= static_cast<int32_t>(field.cols) ||
                pos.y() >= static_cast<int32_t>(field.rows) ||
                pos.x() < 0 || pos.y() < 0) {
                return false;
            }

            // Create a unique hash for the position
            uint64_t hash = (static_cast<uint64_t>(pos.x()) << 32) | static_cast<uint64_t>(pos.y());

            // Check if this position is already used
            if (positionSet.find(hash) != positionSet.end()) {
                return false;
            }

            positionSet.insert(hash);
        }

        return true;
    }

    // Generate a 2D representation of the placement
    Matrix<CellID> getPlacementMatrix(const Scheme& scheme) const {
        const auto& field = scheme.getField();
        const auto& cells = scheme.getCells();

        // Create a matrix filled with invalid cells
        Matrix<CellID> matrix(field.rows, field.cols, Cell::INVALID_ID);

        // Place cells
        for (size_t i = 0; i < cells.size(); ++i) {
            const Coordinate& pos = m_positions[i];
            if (pos.isValid()) {
                matrix[pos.y()][pos.x()] = cells[i].getID();
            }
        }

        // Place fillers
        for (const auto& pos : m_fillers) {
            matrix[pos.y()][pos.x()] = Cell::FILLER_ID;
        }

        return matrix;
    }

    // Generate a 2D representation with actual cell objects
    Matrix<Cell> getCellMatrix(const Scheme& scheme) const {
        const auto& field = scheme.getField();
        const auto& cells = scheme.getCells();

        // Create a matrix filled with invalid cells
        Matrix<Cell> matrix(field.rows, field.cols, Cell(Cell::INVALID_ID));

        // Place cells
        for (size_t i = 0; i < cells.size(); ++i) {
            const Coordinate& pos = m_positions[i];
            if (pos.isValid()) {
                matrix[pos.y()][pos.x()] = cells[i];
            }
        }

        // Place fillers
        Cell fillerCell(Cell::FILLER_ID, "Filler");
        for (const auto& pos : m_fillers) {
            matrix[pos.y()][pos.x()] = fillerCell;
        }

        return matrix;
    }

    // Convert to JSON
    nlohmann::json toJson() const {
        nlohmann::json j;

        // Add cell positions
        nlohmann::json positionsJson = nlohmann::json::array();
        for (const auto& pos : m_positions) {
            positionsJson.push_back(pos.toJson());
        }
        j["positions"] = positionsJson;

        // Add filler positions
        nlohmann::json fillersJson = nlohmann::json::array();
        for (const auto& pos : m_fillers) {
            fillersJson.push_back(pos.toJson());
        }
        j["fillers"] = fillersJson;

        return j;
    }

    // Create from JSON
    static Chromosome fromJson(const nlohmann::json& j, const Scheme& scheme) {
        Chromosome chromosome(scheme);

        // Load cell positions
        if (j.contains("positions") && j["positions"].is_array()) {
            const auto& positionsJson = j["positions"];
            for (size_t i = 0; i < std::min(chromosome.size(), positionsJson.size()); ++i) {
                chromosome[i] = Coordinate::fromJson(positionsJson[i]);
            }
        }

        // Load filler positions
        if (j.contains("fillers") && j["fillers"].is_array()) {
            const auto& fillersJson = j["fillers"];
            chromosome.m_fillers.clear();
            for (const auto& fillerJson : fillersJson) {
                chromosome.m_fillers.push_back(Coordinate::fromJson(fillerJson));
            }
        }

        return chromosome;
    }
};