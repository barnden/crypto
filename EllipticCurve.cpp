#include "EllipticCurve.h"
#include "Modmath.h"

#include <iostream>

void EllipticCurve::GeneratePoints()
{
    for (auto y = 0; y < m_field; y++) {
        // Find quadratic residue in Fm
        auto y2 = Modexp(y, 2, m_field);

        for (auto x = 0; x < m_field; x++) {
            auto value = (x * (x * x + m_a) + m_b) % m_field;

            if (y2 != value)
                continue;

            m_points.push_back(Point(x, y2, m_a, m_b, m_field));
        }
    }
}

bool Point::IsOnCurve()
{
    auto y2 = Modexp(GetY(), 2, GetField());
    auto val = GetX() * (Modexp(GetX(), 2, GetField()) + GetA()) + GetB();

    val %= GetField();

    return y2 == val;
}

void Point::CheckValidity()
{
    assert(IsOnCurve());
}

Point operator-(Point const& rhs)
{
    auto result = rhs;

    result.m_coord.m_y = rhs.GetField() - rhs.GetY();

    return result;
}

Point& Point::operator+=(Point const& rhs)
{
    if (*this ^= rhs) {
        *this = Point(Coordinate { 0, 0, 2 }, 0, 0, 0);
        return *this;
    }

    // If both points are inf, return inf
    if (GetW() == 0 && rhs.GetW() == 0)
        return Point::set_point_at_infinity(*this);

    // If either are points at infinity, then return the other
    if (GetW() == 0) {
        *this = rhs;
        return *this;
    }

    if (rhs.GetW() == 0)
        return *this;

    // Given points L, R on E(F_p), if L.x == R.x return inf
    if (*this != rhs && GetX() == rhs.GetX())
        return Point::set_point_at_infinity(*this);

    auto lambda = 0ull;
    auto const field = GetField();

    auto const lx = GetX();
    auto const ly = GetY();

    auto const rx = rhs.GetX();
    auto const ry = rhs.GetY();

    if (*this == rhs) // Point doubling
        lambda = (3 * (lx * lx) + GetA()) * Modinv(2 * ly, field);
    else
        lambda = Modsub(ry, ly, field) * Modinv(Modsub(rx, lx, field), field);

    lambda %= field;

    auto xn = Modexp(lambda, 2, field);
    xn = Modsub(xn, lx + rx, field);

    auto yn = lambda * Modsub(lx, xn, field);
    yn = Modsub(yn, ly, field);

    // Check if on curve

    m_coord = Coordinate { xn, yn };

    if (IsOnCurve())
        return *this;

    return Point::set_point_at_infinity(*this);
}

Point& Point::operator-=(Point const& rhs)
{
    if (rhs.GetW() == 0)
        return *this;

    return *this += -rhs;
}

Point Point::operator+(Point const& rhs)
{
    auto result = *this;

    return result += rhs;
}

Point Point::operator-(Point const& rhs)
{
    auto result = *this;

    return result -= rhs;
}

Point& Point::operator*=(int64_t rhs)
{
    if (rhs == 1)
        return *this;

    if (rhs == -1) {
        *this = -*this;
        return *this;
    }

    auto nrhs = rhs;
    auto P = *this;

    if (rhs < 0) {
        nrhs *= -1;
        P = -*this;
    }

    auto const m = (sizeof(rhs) * 8) - __builtin_clzll(nrhs);
    Point::set_point_at_infinity(*this);

#ifdef MONTGOMERY
    // Perform the Montgomery ladder method of point addition
    // Computes point multiplication in fixed time
    for (auto i = m; i-- > 0;) {
        if (nrhs & (1 << i)) {
            *this += P;
            P += P;

            continue;
        }

        P += *this;
        *this += *this;
    }
#else
    // Similar to fast powering impl in Modexp from Modmath.h
    // Vulnerable to side-channel attacks
    if (rhs == 0)
        return *this;

    for (auto i = m; i-- > 0;) {
        *this += *this;

        if (nrhs & (1 << i))
            *this += P;
    }
#endif

    return *this;
}

Point operator*(int64_t lhs, Point const& rhs)
{
    auto result = rhs;

    return result *= lhs;
}

Point operator*(Point const& lhs, int rhs) { return rhs * lhs; }

// return true if points are on same curve
inline bool Point::operator|=(Point const& rhs) { return static_cast<Curve>(*this) == static_cast<Curve>(rhs); }
inline bool Point::operator^=(Point const& rhs) { return !(*this |= rhs); }

inline bool Point::operator==(Point const& rhs)
{
    if (*this ^= rhs)
        return false;

    return (GetX() == rhs.GetX()) && (GetY() == rhs.GetY());
}
inline bool Point::operator!=(Point const& rhs) { return !(*this == rhs); }

std::ostream& operator<<(std::ostream& stream, Point const& point)
{
    if (point.GetW() == 0)
        return stream << "inf";

    if (point.GetW() == 2)
        return stream << "err";

    return stream << "(" << point.GetX() << ", " << point.GetY() << ")";
}
