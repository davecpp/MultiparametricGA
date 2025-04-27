#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <memory>
#include <nlohmann/json.hpp>

// Type definitions for clarity
using CellID = int32_t;
using Point = std::pair<int32_t, int32_t>;
using Polygon = std::vector<Point>;
using ConnectionWeight = double;

// Class representing a cell in the IC layout
class Cell {
public:
	// Constants for special cell types
	static constexpr CellID FILLER_ID = -1;
	static constexpr CellID INVALID_ID = -2;

private:
	CellID m_id = INVALID_ID;
	std::string m_name;
	Polygon m_polygon;  // Relative to cell origin
	double m_thermalValue = 0.0;
	double m_powerDensity = 0.0;

public:
	Cell() = default;

	Cell(CellID id, const std::string& name = "")
		: m_id(id), m_name(name) {
	}

	Cell(CellID id, const std::string& name, const Polygon& polygon)
		: m_id(id), m_name(name), m_polygon(polygon) {
	}

	// Getters
	CellID getID() const { return m_id; }
	const std::string& getName() const { return m_name; }
	const Polygon& getPolygon() const { return m_polygon; }
	double getThermalValue() const { return m_thermalValue; }
	double getPowerDensity() const { return m_powerDensity; }

	// Setters
	void setID(CellID id) { m_id = id; }
	void setName(const std::string& name) { m_name = name; }
	void setPolygon(const Polygon& polygon) { m_polygon = polygon; }
	void setThermalValue(double value) { m_thermalValue = value; }
	void setPowerDensity(double value) { m_powerDensity = value; }

	// Helper methods
	bool isValid() const { return m_id != INVALID_ID; }
	bool isFiller() const { return m_id == FILLER_ID; }

	// Get absolute polygon coordinates based on cell position
	Polygon getAbsolutePolygon(int32_t x, int32_t y) const {
		Polygon absolutePolygon;
		absolutePolygon.reserve(m_polygon.size());

		for (const auto& point : m_polygon) {
			absolutePolygon.emplace_back(point.first + x, point.second + y);
		}

		return absolutePolygon;
	}

	// Get bounding box (min_x, min_y, max_x, max_y)
	std::tuple<int32_t, int32_t, int32_t, int32_t> getBoundingBox() const {
		if (m_polygon.empty()) {
			return { 0, 0, 0, 0 };
		}

		int32_t min_x = m_polygon[0].first;
		int32_t min_y = m_polygon[0].second;
		int32_t max_x = m_polygon[0].first;
		int32_t max_y = m_polygon[0].second;

		for (const auto& point : m_polygon) {
			min_x = std::min(min_x, point.first);
			min_y = std::min(min_y, point.second);
			max_x = std::max(max_x, point.first);
			max_y = std::max(max_y, point.second);
		}

		return { min_x, min_y, max_x, max_y };
	}

	// Get cell width and height based on bounding box
	std::pair<int32_t, int32_t> getSize() const {
		auto [min_x, min_y, max_x, max_y] = getBoundingBox();
		return { max_x - min_x, max_y - min_y };
	}

	// Serialize to JSON
	nlohmann::json toJson() const {
		nlohmann::json j;
		j["id"] = m_id;
		j["name"] = m_name;

		// Convert polygon to JSON array
		nlohmann::json polygon = nlohmann::json::array();
		for (const auto& point : m_polygon) {
			polygon.push_back({ point.first, point.second });
		}
		j["polygon"] = polygon;

		j["thermal_value"] = m_thermalValue;
		j["power_density"] = m_powerDensity;

		return j;
	}

	// Deserialize from JSON
	static Cell fromJson(const nlohmann::json& j) {
		Cell cell;

		if (j.contains("id")) {
			cell.m_id = j["id"].get<CellID>();
		}

		if (j.contains("name")) {
			cell.m_name = j["name"].get<std::string>();
		}

		if (j.contains("polygon") && j["polygon"].is_array()) {
			Polygon polygon;
			for (const auto& point : j["polygon"]) {
				if (point.is_array() && point.size() == 2) {
					polygon.emplace_back(point[0].get<int32_t>(), point[1].get<int32_t>());
				}
			}
			cell.m_polygon = polygon;
		}

		if (j.contains("thermal_value")) {
			cell.m_thermalValue = j["thermal_value"].get<double>();
		}

		if (j.contains("power_density")) {
			cell.m_powerDensity = j["power_density"].get<double>();
		}

		return cell;
	}
};

// Type alias for a collection of cells
using CellsVector = std::vector<Cell>;

// 2D template class for matrices
template<typename T>
class Matrix {
private:
	std::vector<std::vector<T>> m_data;
	size_t m_rows = 0;
	size_t m_cols = 0;

public:
	Matrix() = default;

	Matrix(size_t rows, size_t cols)
		: m_rows(rows), m_cols(cols) {
		m_data.resize(rows, std::vector<T>(cols));
	}

	Matrix(size_t rows, size_t cols, const T& defaultValue)
		: m_rows(rows), m_cols(cols) {
		m_data.resize(rows, std::vector<T>(cols, defaultValue));
	}

	// Access operators
	std::vector<T>& operator[](size_t row) {
		return m_data[row];
	}

	const std::vector<T>& operator[](size_t row) const {
		return m_data[row];
	}

	// Size information
	size_t rows() const { return m_rows; }
	size_t cols() const { return m_cols; }

	// Resize the matrix
	void resize(size_t rows, size_t cols) {
		m_rows = rows;
		m_cols = cols;
		m_data.resize(rows);
		for (auto& row : m_data) {
			row.resize(cols);
		}
	}

	// Resize with default value
	void resize(size_t rows, size_t cols, const T& defaultValue) {
		m_rows = rows;
		m_cols = cols;
		m_data.resize(rows);
		for (auto& row : m_data) {
			row.resize(cols, defaultValue);
		}
	}

	// Fill the matrix with a value
	void fill(const T& value) {
		for (auto& row : m_data) {
			std::fill(row.begin(), row.end(), value);
		}
	}

	// Convert to JSON
	nlohmann::json toJson() const {
		return m_data;
	}

	// Create from JSON
	static Matrix<T> fromJson(const nlohmann::json& j) {
		if (!j.is_array()) {
			return Matrix<T>();
		}

		size_t rows = j.size();
		size_t cols = 0;
		if (rows > 0 && j[0].is_array()) {
			cols = j[0].size();
		}

		Matrix<T> matrix(rows, cols);

		for (size_t i = 0; i < rows; ++i) {
			if (j[i].is_array()) {
				for (size_t k = 0; k < std::min(cols, j[i].size()); ++k) {
					matrix[i][k] = j[i][k].get<T>();
				}
			}
		}

		return matrix;
	}
};

// Type alias for connection matrix
using ConnectionMatrix = Matrix<ConnectionWeight>;