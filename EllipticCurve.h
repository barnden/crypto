#pragma once

#define ASSERT_NOT_REACHED __builtin_unreachable()

#include <cstdint>
#include <ostream>
#include <vector>

struct Coordinate {
    // Homogeneous coordinate of a point
    uint64_t m_x;
    uint64_t m_y;
    uint64_t m_w = 1;
};

class EllipticCurve;

class Curve {
public:
    Curve() { ASSERT_NOT_REACHED; }
    Curve(uint64_t a, uint64_t b, uint64_t field)
        : m_field(field)
        , m_a(a)
        , m_b(b)
    {
    }

    inline uint64_t GetField() const { return m_field; };
    inline uint64_t GetA() const { return m_a; };
    inline uint64_t GetB() const { return m_b; };

    inline friend bool operator==(Curve const& lhs, Curve const& rhs)
    {
        return lhs.m_field == rhs.m_field
               && lhs.m_a == rhs.m_a
               && lhs.m_b == rhs.m_b;
    }

    inline friend bool operator!=(Curve const& lhs, Curve const& rhs) { return !(lhs == rhs); }

private:
    friend class EllipticCurve;
    uint64_t m_field;
    uint64_t m_a;
    uint64_t m_b;
};

class Point : public Curve {
public:
    Point(uint64_t x, uint64_t y, uint64_t a, uint64_t b, uint64_t field)
        : m_coord(Coordinate(x, y))
        , Curve(a, b, field)
    {
        CheckValidity();
    }

    Point(uint64_t x, uint64_t y, Curve const& curve)
        : m_coord(Coordinate(x, y))
        , Curve(curve)
    {
        CheckValidity();
    }

    Point(Coordinate coord, uint64_t a, uint64_t b, uint64_t field)
        : m_coord(coord)
        , Curve(a, b, field)
    {
        CheckValidity();
    }

    Point(Coordinate coord, Curve const& curve)
        : m_coord(coord)
        , Curve(curve)
    {
        CheckValidity();
    }

    template <typename T>
    static inline Point make_point_at_infinity(T curve)
    {
        return Point { Coordinate { 0, 0, 0 }, static_cast<Curve>(curve) };
    }

    uint64_t GetX() const { return m_coord.m_x; }
    uint64_t GetY() const { return m_coord.m_y; }

    // The point at infinity is represented with homogeneous coordinate of w=0
    // Points on the curve have w=1
    // For purposes of errors we have w=2
    uint64_t GetW() const { return m_coord.m_w; }
    Coordinate GetCoordinate() const { return m_coord; }

    inline friend bool operator|=(Point const& lhs, Point const& rhs)
    {
        // return true if points are on same curve
        return static_cast<Curve>(lhs) == static_cast<Curve>(rhs);
    }

    inline friend bool operator^=(Point const& lhs, Point const& rhs) { return !(lhs |= rhs); }

    inline friend bool operator==(Point const& lhs, Point const& rhs)
    {
        if (lhs ^= rhs)
            return false;

        return (lhs.GetX() == rhs.GetX()) && (lhs.GetY() == rhs.GetY());
    }

    inline friend bool operator!=(Point const& lhs, Point const& rhs) { return !(lhs == rhs); }

    inline friend std::ostream& operator<<(std::ostream& stream, Point const& point)
    {
        if (point.GetW() == 0)
            return stream << "inf";

        if (point.GetW() == 2)
            return stream << "err";

        return stream << "(" << point.GetX() << ", " << point.GetY() << ")";
    }

private:
    Coordinate m_coord;

    void CheckValidity();
};

class EllipticCurve : public Curve {
public:
    EllipticCurve() { ASSERT_NOT_REACHED; }

    EllipticCurve(uint64_t a, uint64_t b, uint64_t field)
        : Curve(a, b, field)
    {
    }

    inline std::vector<Point> GetPoints() {
        if (!m_points.size())
            GeneratePoints();

        return m_points;
    };

private:
    std::vector<Point> m_points;

    void GeneratePoints();
};

Point operator+(Point lhs, Point const& rhs);
Point operator-(Point lhs, Point const& rhs);
Point operator-(Point lhs);
Point operator*(int64_t lhs, Point const& rhs);
Point operator*(Point const& lhs, int64_t rhs);
Point& operator+=(Point& lhs, Point const& rhs);
