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

void Point::CheckValidity()
{
    auto y2 = Modexp(GetY(), 2, GetField());
    auto x = GetX();
    auto value = (x * (x * x + GetA()) + GetB()) % GetField();

    if (y2 != value)
        ASSERT_NOT_REACHED;
}

Point operator+(Point lhs, Point const& rhs)
{
    if (lhs ^= rhs)
        return Point(Coordinate { 0, 0, 2 }, 0, 0, 0);

    // If both are points at infinity, return infty
    // If L.x == R.x, return infty
    if (lhs.GetW() == 0 && rhs.GetW() == 0)
        return Point::make_point_at_infinity(lhs);

    // If either are points at infinity, then return the other
    if (lhs.GetW() == 0)
        return rhs;

    if (rhs.GetW() == 0)
        return lhs;

    // Given L, R in E(F_p), if L != R and L.x == R.x, then return infty
    if (lhs != rhs && lhs.GetX() == rhs.GetX())
        return Point::make_point_at_infinity(lhs);

    auto lambda = 0ull;
    auto const field = lhs.GetField();

    auto const lx = lhs.GetX();
    auto const ly = lhs.GetY();

    auto const rx = rhs.GetX();
    auto const ry = rhs.GetY();

    if (lhs == rhs) // Point doubling
        lambda = (3 * (lx * lx) + lhs.GetA()) * Modinv(2 * ly, field);
    else
        lambda = Modsub(ry, ly, field) * Modinv(Modsub(rx, lx, field), field);

    lambda %= field;

    auto xn = Modexp(lambda, 2, field);
    xn = Modsub(xn, lx + rx, field);

    auto yn = lambda * Modsub(lx, xn, field);
    yn = Modsub(yn, ly, field);

    xn %= field;
    yn %= field;

    // Check if on curve
    auto y2 = Modexp(yn, 2, field);
    auto val = xn * (Modexp(xn, 2, field) + lhs.GetA()) + lhs.GetB();

    val %= field;

    if (y2 != val)
        return Point::make_point_at_infinity(lhs);

    return Point(Coordinate { xn, yn }, static_cast<Curve>(lhs));
}

Point operator-(Point lhs, Point const& rhs)
{
    if (rhs.GetW() == 0)
        return lhs;

    return lhs + (-rhs);
}

Point operator-(Point lhs)
{
    return Point {
        Coordinate { lhs.GetX(), lhs.GetField() - lhs.GetY(), lhs.GetW() },
        static_cast<Curve>(lhs)
    };
}

Point operator*(int64_t lhs, Point const& rhs)
{
    if (lhs == 1)
        return rhs;

    if (lhs == -1)
        return -rhs;

    auto nlhs = lhs;
    auto P = rhs;

    if (lhs < 0) {
        nlhs = -nlhs;
        P = -rhs;
    }

    auto const m = (sizeof(lhs) * 8) - __builtin_clzll(nlhs);
    auto accumulator = Point::make_point_at_infinity(rhs);

#ifdef MONTGOMERY
    // Perform the Montgomery ladder method of point addition
    // Computes point multiplication in fixed time
    for (auto i = m; i-- > 0;) {
        if (nlhs & (1 << i)) {
            accumulator += P;
            P += P;

            continue;
        }

        P += accumulator;
        accumulator += accumulator;
    }

    return accumulator;
#else
    // Similar to fast powering impl in Modexp from Modmath.h
    // Vulnerable to side-channel attacks
    if (lhs == 0)
        return accumulator;

    for (auto i = m; i-- > 0;) {
        accumulator += accumulator;

        if (nlhs & (1 << i))
            accumulator += P;
    }

    return accumulator;
#endif
}

Point operator*(Point const& lhs, int rhs) { return rhs * lhs; }

Point& operator+=(Point& lhs, Point const& rhs)
{
    lhs = lhs + rhs;

    return lhs;
}
