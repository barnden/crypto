#pragma once

#include <cassert>
#include <cstdint>
#include <ostream>
#include <vector>

#include "BigInt.h"

#define ASSERT_NOT_REACHED assert(false)

struct Coordinate {
    // Homogeneous coordinate of a point
    BigInt m_x;
    BigInt m_y;
    bool m_w = 1; // w is 0 for point at infinity; 1 otherwise
};

class EllipticCurve;

class Curve {
public:
    Curve() { ASSERT_NOT_REACHED; }
    Curve(BigInt a, BigInt b, BigInt field)
        : m_field(field)
        , m_a(a)
        , m_b(b)
    {
    }

    inline BigInt get_field() const { return m_field; };
    inline BigInt get_a() const { return m_a; };
    inline BigInt get_b() const { return m_b; };

    inline friend bool operator==(Curve const& lhs, Curve const& rhs)
    {
        return lhs.m_field == rhs.m_field
               && lhs.m_a == rhs.m_a
               && lhs.m_b == rhs.m_b;
    }

    inline friend bool operator!=(Curve const& lhs, Curve const& rhs) { return !(lhs == rhs); }

private:
    friend class EllipticCurve;
    BigInt m_field;
    BigInt m_a;
    BigInt m_b;
};

class Point : public Curve {
public:
    Point(BigInt x, BigInt y, BigInt a, BigInt b, BigInt field)
        : m_coord(Coordinate(x, y))
        , Curve(a, b, field)
    {
        check_validity();
    }

    Point(BigInt x, BigInt y, Curve const& curve)
        : m_coord(Coordinate(x, y))
        , Curve(curve)
    {
        check_validity();
    }

    Point(Coordinate coord, BigInt a, BigInt b, BigInt field)
        : m_coord(coord)
        , Curve(a, b, field)
    {
        check_validity();
    }

    Point(Coordinate coord, Curve const& curve)
        : m_coord(coord)
        , Curve(curve)
    {
        check_validity();
    }

    template <typename T>
    static inline Point make_point_at_infinity(T curve)
    {
        return Point { Coordinate { 0, 0, 0 }, static_cast<Curve>(curve) };
    }

    static inline Point& set_point_at_infinity(Point& point)
    {
        point.m_coord = Coordinate { 0, 0, 0 };
        return point;
    }

    BigInt get_x() const { return m_coord.m_x; }
    BigInt get_y() const { return m_coord.m_y; }

    // The point at infinity is represented with homogeneous coordinate of w=0; and w=1 otherwise
    bool get_w() const { return m_coord.m_w; }
    Coordinate get_coordinate() const { return m_coord; }

    Point& operator+=(Point const& rhs);
    Point& operator-=(Point const& rhs);
    Point& operator*=(int64_t rhs);

    Point operator+(Point const& rhs) const;
    Point operator-(Point const& rhs) const;
    Point operator-() const;

    friend Point operator*(int64_t lhs, Point const& rhs);
    friend Point operator*(Point const& lhs, int64_t rhs);

    bool operator|=(Point const& rhs);
    bool operator^=(Point const& rhs);
    bool operator==(Point const& rhs);
    bool operator!=(Point const& rhs);

    friend std::ostream& operator<<(std::ostream& stream, Point const& point);

private:
    Coordinate m_coord;

    bool is_on_curve();
    void check_validity();
};

class EllipticCurve : public Curve {
private:
    std::vector<Point> m_points;

    void generate_points();

public:
    EllipticCurve() { ASSERT_NOT_REACHED; }

    EllipticCurve(BigInt a, BigInt b, BigInt field)
        : Curve(a, b, field)
    {
    }

    inline std::vector<Point> get_points()
    {
        if (!m_points.size())
            generate_points();

        return m_points;
    };
};
