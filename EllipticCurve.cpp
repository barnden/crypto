#include "EllipticCurve.h"
#include "Modmath.h"

void EllipticCurve::generate_points()
{
    for (auto y = BigInt {}; y < m_field; y += 1) {
        // Find quadratic residues in m_field
        auto y2 = Modexp(y, 2, m_field);

        for (auto x = BigInt {}; x < m_field; x += 1) {
            auto value = (x * (x * x + m_a) + m_b) % m_field;

            if (y2 == value)
                m_points.push_back(Point(x, y, m_a, m_b, m_field));
        }
    }
}

bool Point::is_on_curve()
{
    if (get_w() == 0)
        return true;

    auto y2 = Modexp(get_y(), 2, get_field());
    auto val = (get_x() * (get_x() * get_x() + get_a()) + get_b()) % get_field();

    return y2 == val;
}

void Point::check_validity()
{
    assert(is_on_curve());
}

Point Point::operator-() const
{
    auto result = Point { *this };

    result.m_coord.m_y = get_field() - get_y();

    return result;
}

Point& Point::operator+=(Point const& rhs)
{
    assert(*this |= rhs);

    // If both points are inf, return inf
    if (get_w() == 0 && rhs.get_w() == 0)
        return Point::set_point_at_infinity(*this);

    // If either are points at infinity, then return the other
    if (get_w() == 0) {
        *this = rhs;
        return *this;
    }

    if (rhs.get_w() == 0)
        return *this;

    // Given points L, R on E(F_p), if L.x == R.x return inf
    if (*this != rhs && get_x() == rhs.get_x())
        return Point::set_point_at_infinity(*this);

    auto lambda = BigInt {};
    auto const field = get_field();

    auto const lx = get_x();
    auto const ly = get_y();

    auto const rx = rhs.get_x();
    auto const ry = rhs.get_y();

    if (*this == rhs) // Point doubling
        lambda = ((lx * lx) * 3 + get_a()) * Modinv(ly * 2, field);
    else
        lambda = ((ry - ly) % field) * Modinv((rx - lx) % field, field);

    lambda %= field;

    auto xn = (Modexp(lambda, 2, field) - (lx + rx)) % field;
    auto yn = (lambda * (lx - xn) - ly) % field;

    // Check if on curve

    m_coord = Coordinate { xn, yn };

    if (is_on_curve())
        return *this;

    return Point::set_point_at_infinity(*this);
}

Point& Point::operator-=(Point const& rhs)
{
    if (rhs.get_w() == 0)
        return *this;

    return *this += -rhs;
}

Point Point::operator+(Point const& rhs) const
{
    auto result = Point { *this };

    return result += rhs;
}

Point Point::operator-(Point const& rhs) const
{
    auto result = Point { *this };

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

Point operator*(Point const& lhs, int64_t rhs) { return rhs * lhs; }

// return true if points are on same curve
inline bool Point::operator|=(Point const& rhs) { return static_cast<Curve>(*this) == static_cast<Curve>(rhs); }
inline bool Point::operator^=(Point const& rhs) { return !(*this |= rhs); }

inline bool Point::operator==(Point const& rhs)
{
    if (*this ^= rhs)
        return false;

    if (get_w() == 0 && rhs.get_w() == 0)
        return true;

    return (get_x() == rhs.get_x()) && (get_y() == rhs.get_y());
}
inline bool Point::operator!=(Point const& rhs) { return !(*this == rhs); }

std::ostream& operator<<(std::ostream& stream, Point const& point)
{
    if (point.get_w() == 0)
        return stream << "inf";

    return stream << "(" << point.get_x() << ", " << point.get_y() << ")";
}
