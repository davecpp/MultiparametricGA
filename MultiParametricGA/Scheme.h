#pragma once

#include "Cell.h"
#include <fstream>
#include <iostream>
#include <string>
#include <nlohmann/json.hpp>

// Class representing a placement field configuration
struct PlacementField {
    size_t rows = 0;
    size_t cols = 0;
    bool allowFillers = true;

    // Convert to JSON
    nlohmann::json toJson() const {
        nlohmann::json j;
        j["rows"] = rows;
        j["cols"] = cols;
        j["allow_fillers"] = allowFillers;
        return j;
    }

    // Create from JSON
    static PlacementField fromJson(const nlohmann::json& j) {
        PlacementField field;

        if (j.contains("rows")) {
            field.rows = j["rows"].get<size_t>();
        }

        if (j.contains("cols")) {
            field.cols = j["cols"].get<size_t>();
        }

        if (j.contains("allow_fillers")) {
            field.allowFillers = j["allow_fillers"].get<bool>();
        }

        return field;
    }
};

// Class representing a complete IC scheme with cells and connections
class Scheme {
private:
    CellsVector m_cells;
    ConnectionMatrix m_connections;
    PlacementField m_field;

public:
    Scheme() = default;

    // Getters
    const CellsVector& getCells() const { return m_cells; }
    const ConnectionMatrix& getConnections() const { return m_connections; }
    const PlacementField& getField() const { return m_field; }

    // Setters
    void setCells(const CellsVector& cells) { m_cells = cells; }
    void setConnections(const ConnectionMatrix& connections) { m_connections = connections; }
    void setField(const PlacementField& field) { m_field = field; }

    // Add a single cell
    void addCell(const Cell& cell) {
        m_cells.push_back(cell);
    }

    // Get cell by ID
    Cell* getCellById(CellID id) {
        for (auto& cell : m_cells) {
            if (cell.getID() == id) {
                return &cell;
            }
        }
        return nullptr;
    }

    // Check if the scheme is valid
    bool isValid() const {
        // Check if field is properly sized
        if (m_field.rows == 0 || m_field.cols == 0) {
            return false;
        }

        // Check if there are any cells
        if (m_cells.empty()) {
            return false;
        }

        // Check if connection matrix matches the number of cells
        if (m_connections.rows() != m_cells.size() || m_connections.cols() != m_cells.size()) {
            return false;
        }

        return true;
    }

    // Load scheme from JSON file
    bool loadFromFile(const std::string& filename) {
        try {
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Error: Unable to open scheme file: " << filename << std::endl;
                return false;
            }

            nlohmann::json j;
            file >> j;

            return loadFromJson(j);
        }
        catch (const std::exception& e) {
            std::cerr << "Error loading scheme from file: " << e.what() << std::endl;
            return false;
        }
    }

    // Save scheme to JSON file
    bool saveToFile(const std::string& filename) const {
        try {
            std::ofstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
                return false;
            }

            nlohmann::json j = toJson();
            file << j.dump(4); // Pretty print with 4 spaces
            return true;
        }
        catch (const std::exception& e) {
            std::cerr << "Error saving scheme to file: " << e.what() << std::endl;
            return false;
        }
    }

    // Convert to JSON
    nlohmann::json toJson() const {
        nlohmann::json j;

        // Add cells
        nlohmann::json cellsJson = nlohmann::json::array();
        for (const auto& cell : m_cells) {
            cellsJson.push_back(cell.toJson());
        }
        j["cells"] = cellsJson;

        // Add connections
        j["connections"] = m_connections.toJson();

        // Add field configuration
        j["field"] = m_field.toJson();

        return j;
    }

    // Load from JSON
    bool loadFromJson(const nlohmann::json& j) {
        try {
            // Load cells
            if (j.contains("cells") && j["cells"].is_array()) {
                m_cells.clear();
                for (const auto& cellJson : j["cells"]) {
                    m_cells.push_back(Cell::fromJson(cellJson));
                }
            }

            // Load connections
            if (j.contains("connections") && j["connections"].is_array()) {
                m_connections = ConnectionMatrix::fromJson(j["connections"]);
            }

            // Load field configuration
            if (j.contains("field") && j["field"].is_object()) {
                m_field = PlacementField::fromJson(j["field"]);
            }

            // Ensure connections matrix has correct size
            if (m_connections.rows() != m_cells.size() || m_connections.cols() != m_cells.size()) {
                m_connections.resize(m_cells.size(), m_cells.size(), 0.0);
            }

            return true;
        }
        catch (const std::exception& e) {
            std::cerr << "Error parsing scheme JSON: " << e.what() << std::endl;
            return false;
        }
    }

    // Generate a test scheme with random cells and connections
    static Scheme generateTestScheme(size_t cellCount, size_t fieldRows, size_t fieldCols) {
        Scheme scheme;

        // Set field
        PlacementField field;
        field.rows = fieldRows;
        field.cols = fieldCols;
        field.allowFillers = true;
        scheme.setField(field);

        // Generate cells
        CellsVector cells;
        for (size_t i = 0; i < cellCount; ++i) {
            Cell cell(static_cast<CellID>(i), "Cell_" + std::to_string(i));

            // Generate a simple rectangular polygon
            int32_t width = 1 + rand() % 3;
            int32_t height = 1 + rand() % 3;

            Polygon polygon = {
                {0, 0},
                {width, 0},
                {width, height},
                {0, height},
                {0, 0}  // Close the polygon
            };

            cell.setPolygon(polygon);
            cell.setThermalValue(0.1 + (rand() % 100) / 100.0);
            cell.setPowerDensity(0.05 + (rand() % 100) / 200.0);

            cells.push_back(cell);
        }
        scheme.setCells(cells);

        // Generate connection matrix
        ConnectionMatrix connections(cellCount, cellCount, 0.0);
        for (size_t i = 0; i < cellCount; ++i) {
            for (size_t j = i + 1; j < cellCount; ++j) {
                // Generate a random connection weight
                double weight = (rand() % 100) / 100.0;
                connections[i][j] = weight;
                connections[j][i] = weight; // Make it symmetric
            }
        }
        scheme.setConnections(connections);

        return scheme;
    }
};