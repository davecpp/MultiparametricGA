#pragma once

#include "Scheme.h"
#include "Chromosome.h"
#include "Cell.h"
#include <string>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

// Class for exporting placement results to GeoJSON
class GeoJSONExporter {
public:
	// Export a placed scheme to GeoJSON
	static nlohmann::json exportToGeoJSON(const Chromosome& chromosome, const Scheme& scheme) {
		nlohmann::json geojson;
		geojson["type"] = "FeatureCollection";
		geojson["features"] = nlohmann::json::array();

		// Get the cell matrix
		Matrix<Cell> cellMatrix = chromosome.getCellMatrix(scheme);
		const auto& field = scheme.getField();
		const auto& cells = scheme.getCells();

		// Add each cell as a feature
		for (size_t y = 0; y < field.rows; ++y) {
			for (size_t x = 0; x < field.cols; ++x) {
				const Cell& cell = cellMatrix[y][x];

				// Skip invalid cells
				if (cell.getID() == Cell::INVALID_ID) {
					continue;
				}

				// Add cell as a feature
				geojson["features"].push_back(cellToFeature(cell, x, y));
			}
		}

		// Add connections as lines
		const auto& connections = scheme.getConnections();
		for (size_t i = 0; i < cells.size(); ++i) {
			for (size_t j = i + 1; j < cells.size(); ++j) {
				double weight = connections[i][j];
				if (weight > 0) {
					const Coordinate& pos1 = chromosome[i];
					const Coordinate& pos2 = chromosome[j];

					// Skip invalid positions
					if (!pos1.isValid() || !pos2.isValid()) {
						continue;
					}

					// Add connection as a line feature
					geojson["features"].push_back(
						connectionToFeature(pos1, pos2, weight, cells[i], cells[j])
					);
				}
			}
		}

		// Add metadata
		geojson["metadata"] = schemeToMetadata(scheme);

		return geojson;
	}

	// Save GeoJSON to file
	static bool saveToFile(const nlohmann::json& geojson, const std::string& filename) {
		try {
			std::ofstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
				return false;
			}

			file << geojson.dump(4); // Pretty print with 4 spaces
			return true;
		}
		catch (const std::exception& e) {
			std::cerr << "Error saving GeoJSON to file: " << e.what() << std::endl;
			return false;
		}
	}

	// Convenience method to export and save in one step
	static bool exportToFile(const Chromosome& chromosome, const Scheme& scheme, const std::string& filename) {
		nlohmann::json geojson = exportToGeoJSON(chromosome, scheme);
		return saveToFile(geojson, filename);
	}

private:
	// Convert a cell to a GeoJSON feature
	static nlohmann::json cellToFeature(const Cell& cell, int32_t x, int32_t y) {
		nlohmann::json feature;
		feature["type"] = "Feature";

		// Properties
		nlohmann::json properties;
		properties["id"] = cell.getID();
		properties["name"] = cell.getName();
		properties["type"] = cell.isFiller() ? "filler" : "cell";
		properties["thermal_value"] = cell.getThermalValue();
		properties["power_density"] = cell.getPowerDensity();
		feature["properties"] = properties;

		// Geometry
		nlohmann::json geometry;
		geometry["type"] = "Polygon";

		// Get absolute polygon coordinates
		Polygon polygon = cell.getAbsolutePolygon(x, y);

		// Convert to GeoJSON format (with correct nesting)
		// Note: GeoJSON format for polygons requires array of arrays of coordinates
		nlohmann::json coordinates = nlohmann::json::array();
		nlohmann::json ring = nlohmann::json::array();

		// Skip empty polygons for fillers
		if (!polygon.empty()) {
			for (const auto& point : polygon) {
				ring.push_back({ point.first, point.second });
			}

			// Ensure polygon is closed
			if (polygon.size() > 0 && (polygon.front().first != polygon.back().first ||
				polygon.front().second != polygon.back().second)) {
				ring.push_back({ polygon.front().first, polygon.front().second });
			}

			coordinates.push_back(ring);
			geometry["coordinates"] = coordinates;
		}
		else {
			// For fillers or empty cells, create an empty array
			geometry["coordinates"] = nlohmann::json::array();
		}

		feature["geometry"] = geometry;

		return feature;
	}

	// Convert a connection to a GeoJSON feature
	static nlohmann::json connectionToFeature(const Coordinate& pos1, const Coordinate& pos2,
		double weight, const Cell& cell1, const Cell& cell2) {
		nlohmann::json feature;
		feature["type"] = "Feature";

		// Properties
		nlohmann::json properties;
		properties["type"] = "connection";
		properties["weight"] = weight;
		properties["source_id"] = cell1.getID();
		properties["target_id"] = cell2.getID();
		properties["source_name"] = cell1.getName();
		properties["target_name"] = cell2.getName();
		feature["properties"] = properties;

		// Geometry
		nlohmann::json geometry;
		geometry["type"] = "LineString";

		// Convert to GeoJSON format
		nlohmann::json coordinates = nlohmann::json::array();
		coordinates.push_back({ pos1.x(), pos1.y() });
		coordinates.push_back({ pos2.x(), pos2.y() });

		geometry["coordinates"] = coordinates;

		feature["geometry"] = geometry;

		return feature;
	}

	// Convert scheme information to metadata
	static nlohmann::json schemeToMetadata(const Scheme& scheme) {
		nlohmann::json metadata;
		metadata["field"] = scheme.getField().toJson();
		metadata["cell_count"] = scheme.getCells().size();

		// Count connections
		const auto& connections = scheme.getConnections();
		size_t connectionCount = 0;
		double totalWeight = 0.0;

		for (size_t i = 0; i < connections.rows(); ++i) {
			for (size_t j = i + 1; j < connections.cols(); ++j) {
				if (connections[i][j] > 0) {
					connectionCount++;
					totalWeight += connections[i][j];
				}
			}
		}

		metadata["connection_count"] = connectionCount;
		metadata["total_connection_weight"] = totalWeight;

		return metadata;
	}
};