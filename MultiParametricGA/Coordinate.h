#pragma once

#include <cstdint>
#include <cmath>
#include <nlohmann/json.hpp>

// Class representing a 2D coordinate
class Coordinate {
private:
    int32_t m_x = -1;
    int32_t m_y = -1;

public:
    // Constants
    static constexpr int32_t INVALID_COORD = -1;

    // Constructors
    Coordinate() = default;

    Coordinate(int32_t x, int32_t y) : m_x(x), m_y(y) {}

    // Getters
    int32_t x() const { return m_x; }
    int32_t y() const { return m_y; }

    // Setters
    void setX(int32_t x) { m_x = x; }
    void setY(int32_t y) { m_y = y; }

    // Check if coordinate is valid
    bool isValid() const {
        return m_x != INVALID_COORD && m_y != INVALID_COORD;
    }

    // Create an invalid coordinate
    static Coordinate invalid() {
        return Coordinate(INVALID_COORD, INVALID_COORD);
    }

    // Equality operators
    bool operator==(const Coordinate& other) const {
        return m_x == other.m_x && m_y == other.m_y;
    }

    bool operator!=(const Coordinate& other) const {
        return !(*this == other);
    }

    // Calculate Manhattan distance between two coordinates
    static int32_t manhattanDistance(const Coordinate& a, const Coordinate& b) {
        return std::abs(a.m_x - b.m_x) + std::abs(a.m_y - b.m_y);
    }

    // Calculate Euclidean distance between two coordinates
    static double euclideanDistance(const Coordinate& a, const Coordinate& b) {
        double dx = static_cast<double>(a.m_x - b.m_x);
        double dy = static_cast<double>(a.m_y - b.m_y);
        return std::sqrt(dx * dx + dy * dy);
    }

    // Default distance calculation method (Manhattan)
    static int32_t distance(const Coordinate& a, const Coordinate& b) {
        return manhattanDistance(a, b);
    }

    // Convert to JSON
    nlohmann::json toJson() const {
        return nlohmann::json::array({ m_x, m_y });
    }

    // Create from JSON
    static Coordinate fromJson(const nlohmann::json& j) {
        if (j.is_array() && j.size() == 2) {
            return Coordinate(j[0].get<int32_t>(), j[1].get<int32_t>());
        }
        return Coordinate::invalid();
    }
};