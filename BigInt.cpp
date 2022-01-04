#include "BigInt.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

void BigInt::embiggen(BigInt const& other)
{
    // Pad m_groups to be other.m_groups.size + 1
    if (groups() > other.groups())
        return;

    embiggen(other.groups() + 1);
}

void BigInt::embiggen(size_t size)
{
    if (groups() > size)
        return;

    for (auto i = 0; i <= size - groups(); i++)
        m_groups.push_back(0);
}

inline void emsmallen(std::deque<uint32_t>& groups)
{
    // Remove all leading digit groups valued at 0, except for the last one
    while (groups.back() == 0 && groups.size() > 1)
        groups.pop_back();
}

[[gnu::flatten]] inline void BigInt::emsmallen() { ::emsmallen(m_groups); }

BigInt naive_multiplication(BigInt const& x, uint64_t y)
{
    auto z = std::deque<uint32_t>(x.groups() + 2);
    auto carry = uint64_t {};

    for (auto i = 0; i < x.groups(); i++) {
        auto product = static_cast<uint64_t>(x.m_groups[i]) * y + z[i];

        z[i] = static_cast<uint32_t>(product);
        z[i + 1] = product >> 32;
    }

    emsmallen(z);

    return { z };
}

BigInt naive_multiplication(BigInt const& x, BigInt const& y)
{
    if (y.groups() == 1)
        return naive_multiplication(x, y.m_groups.back());

    auto z = std::deque<uint32_t>(x.groups() + y.groups() + 1);
    auto carry = uint64_t {};

    for (auto i = 0; i < x.groups(); i++) {
        carry = 0;

        for (auto j = 0, k = i; j < y.groups(); j++, k++) {
            auto product = static_cast<uint64_t>(x.m_groups[i]) * static_cast<uint64_t>(y.m_groups[j]) + z[k] + carry;

            z[k] = static_cast<uint32_t>(product);
            carry = product >> 32;
        }

        z[i + y.groups()] = carry;
    }

    emsmallen(z);

    return { z };
}

BigInt naive_muladd(BigInt const& x, uint64_t mul, uint64_t add)
{
    auto z = naive_multiplication(x, mul);

    auto carry = add;
    for (auto i = 0; i < z.groups(); i++) {
        auto sum = static_cast<uint64_t>(z.m_groups[i]) + carry;
        z.m_groups[i] = static_cast<uint32_t>(sum);
        carry = sum >> 32;
    }

    if (carry)
        z.m_groups.push_back(carry);

    return z;
}

BigInt knuth(BigInt const& x, uint64_t y, bool remainder)
{
    auto Q = std::deque<uint32_t>(x.groups());
    auto k = uint64_t {};

    for (auto j = x.groups(); j-- > 0;) {
        uint64_t t = (k << 32) + x.m_groups[j];
        Q[j] = t / y;
        k = t - Q[j] * y;
    }

    if (remainder)
        return BigInt(std::deque<uint32_t> { static_cast<uint32_t>(k), static_cast<uint32_t>(k >> 32) });

    return BigInt { Q };
}

BigInt knuth(BigInt const& x, BigInt const& y, bool remainder)
{
    // The Art of Computer Programming Vol. 2 Seminumerical Algorithms 3rd ed. pg. 284
    // Hacker's Delight divmnu64.c

    if (y.groups() == 1)
        return knuth(x, y.m_groups[0], remainder);

    auto static constexpr B = 1ull << 32;
    auto Q = std::deque<uint32_t>(x.groups());
    auto S = __builtin_clz(y.m_groups.back());
    auto D = B >> S;
    auto U = naive_multiplication(x, D).m_groups;
    auto V = naive_multiplication(y, D).m_groups;

    emsmallen(U);
    emsmallen(V);

    U.push_back(0); // |U| = m + n + 1

    auto n = V.size();
    auto m = U.size() - n;

    for (auto j = m; j-- > 0;) {
        auto qhat = ((static_cast<uint64_t>(U[n + j]) << 32) ^ U[n + j - 1]) / V[n - 1];
        auto rhat = ((static_cast<uint64_t>(U[n + j]) << 32) ^ U[n + j - 1]) % V[n - 1];

        while (qhat >= B || qhat * V[n - 2] > ((rhat << 32) ^ U[n + j - 2])) {
            qhat--;
            rhat += V[n - 1];

            if (rhat >= B)
                break;
        }

        auto k = uint64_t {};
        auto t = uint64_t {};
        for (auto i = 0; i < n; i++) {
            auto p = qhat * V[i];
            t = U[i + j] - k - static_cast<int64_t>(p);
            U[i + j] = t;
            k = (p >> 32) - (t >> 32);
        }

        t = U[n + j] - k;
        U[n + j] = t;

        Q[j] = qhat;
        if (t < 0) {
            Q[j] = Q[j] - 1;
            k = 0;

            for (auto i = 0; i < n; i++) {
                t = static_cast<uint64_t>(U[i + j]) + V[i] + k;
                U[i + j] = t;
                k = t >> 32;
            }

            U[n + j] += k;
        }
    }

    if (remainder)
        return BigInt { Q } >>= S;

    emsmallen(Q);

    return { Q };
}

BigInt::BigInt()
    : m_negative(false)
{
    m_groups.push_back(0);
}

BigInt::BigInt(uint64_t number)
    : m_negative(false)
{
    auto static constexpr offset = sizeof(uint32_t) * 8;

    if (number < 0) {
        m_negative = true;
        number *= -1;
    }

    while (number) {
        m_groups.push_back(static_cast<uint32_t>(number));

        number >>= offset;
    }
}

BigInt::BigInt(std::string number)
{
    m_negative = number[0] == '-';

    bool skip_first = m_negative || number[0] == '+';

    // Ignore all spaces, commas, single quotes
    [&]<char... I>(std::integer_sequence<char, I...>)
    {
        ([&](char c) {
            std::string::iterator delim;
            while ((delim = std::find(number.begin(), number.end(), c)) != number.end())
                number.erase(delim);
        }(I),
         ...);
    }
    (std::integer_sequence<char, ',', ' ', '\''> {});

    auto size = number.size();
    while (size > skip_first) {
        auto place = 1ull;
        auto num = uint32_t {};
        auto stop = size - digits;

        if (__builtin_sub_overflow_p(size, digits, stop))
            stop = 0;

        for (auto i = size; i-- > stop; size--) {
            auto c = number[i];

            if (c < '0' || c > '9')
                break;

            num += (c - '0') * place;

            place *= 10;
        }

        if (num || !m_groups.size())
            m_groups.push_back(num);
    }

    auto z = BigInt {};

    for (auto git = m_groups.rbegin(); git != m_groups.rend(); git++)
        z = naive_muladd(z, base, *git);

    m_groups = z.m_groups;

    emsmallen();
}

BigInt::BigInt(std::deque<uint32_t> group)
    : m_groups(group)
    , m_negative(false)
{
    emsmallen();
}

size_t BigInt::size() const
{
    // Get size of number in bits
    return ((groups() - 1) * 32) + (32 - __builtin_clz(m_groups.back()));
}

BigInt& BigInt::operator-=(BigInt const& rhs)
{
    if (rhs.m_negative)
        return *this += -rhs;

    embiggen(rhs);

    auto lit = m_groups.begin();
    auto const& lend = m_groups.end();

    auto rit = std::as_const(rhs.m_groups).begin();
    auto const& rend = rhs.m_groups.end();

    auto sum = uint64_t {};
    while (lit != lend || rit != rend) {
        if (lit != lend) {
            sum += *lit;
            lit++;
        }

        if (rit != rend) {
            sum -= *rit;
            rit++;
        }

        if (sum < 0) {
            *std::prev(lit) = -sum;
            sum = 0;
            m_negative = true;
        } else {
            *std::prev(lit) = static_cast<uint32_t>(sum);
            sum >>= 32;

            if (sum != 0)
                m_negative = false;
        }
    }

    emsmallen();

    return *this;
}

BigInt& BigInt::operator+=(BigInt const& rhs)
{
    if (rhs.m_negative)
        return *this -= rhs;

    if (m_negative) {
        auto temp = *this;
        *this = rhs;

        return *this -= -temp;
    }

    embiggen(rhs);

    auto lit = m_groups.begin();
    auto const& lend = m_groups.end();

    auto rit = std::as_const(rhs.m_groups).begin();
    auto const& rend = rhs.m_groups.end();

    auto sum = uint64_t {};
    while (lit != lend || rit != rend) {
        if (lit != lend) {
            sum += *lit;
            lit++;
        }

        if (rit != rend) {
            sum += *rit;
            rit++;
        }

        *std::prev(lit) = static_cast<uint32_t>(sum);
        sum >>= 32;
    }

    if (sum)
        m_groups.push_back(1);

    emsmallen();

    return *this;
}

BigInt& BigInt::operator*=(BigInt const& rhs)
{
    // Perform multiplication
    // TODO: Implement Karatsuba and Toom-k

    m_negative ^= rhs.m_negative;
    m_groups = naive_multiplication(*this, rhs).m_groups;

    emsmallen();

    return *this;
}

template <Numeric T>
BigInt& BigInt::operator*=(T rhs)
{
    if (rhs < 0) {
        m_negative ^= 1;
        rhs *= -1;
    }

    m_groups = naive_multiplication(*this, rhs).m_groups;

    emsmallen();

    return *this;
}

BigInt BigInt::operator*(BigInt const& rhs) const { return BigInt { *this } *= rhs; }

template <Numeric T>
BigInt BigInt::operator*(T rhs) const { return BigInt { *this } *= rhs; }

BigInt& BigInt::operator/=(BigInt const& rhs)
{
    m_negative ^= rhs.m_negative;
    m_groups = knuth(*this, rhs, false).m_groups;

    emsmallen();

    return *this;
}

template <Numeric T>
BigInt& BigInt::operator/=(T rhs)
{
    if (rhs < 0) {
        m_negative ^= 1;
        rhs *= -1;
    }

    m_groups = knuth(*this, rhs, false).m_groups;

    emsmallen();

    return *this;
}

BigInt BigInt::operator/(BigInt const& rhs) const { return BigInt { *this } /= rhs; }

template <Numeric T>
BigInt BigInt::operator/(T rhs) const { return BigInt { *this } /= rhs; }

BigInt& BigInt::operator<<=(int rhs)
{
    if (rhs == 0)
        return *this;

    if (rhs < 0)
        return *this >>= -rhs;

    auto groups = rhs / 32;

    for (auto i = 0; i < groups; i++)
        m_groups.push_front(0);

    auto s = rhs % 32;

    if (__builtin_clz(m_groups.back()) < s)
        m_groups.push_back(0);

    for (auto i = m_groups.size() - 1; i > 0; i--)
        m_groups[i] = (m_groups[i] << s) | (m_groups[i - 1] >> (32 - s));

    m_groups[0] <<= s;

    emsmallen();

    return *this;
}

BigInt& BigInt::operator>>=(int rhs)
{
    if (rhs == 0)
        return *this;

    if (rhs < 0)
        return *this <<= -rhs;

    if (rhs >= size()) {
        m_groups.clear();
        m_groups[0] = 0;

        return *this;
    }

    auto groups = rhs / 32;
    for (auto i = 0; i < groups; i++)
        m_groups.pop_front();

    auto s = rhs % 32;
    for (auto i = 0; i < m_groups.size() - 1; i++)
        m_groups[i] = (m_groups[i] >> s) | static_cast<uint64_t>(m_groups[i + 1] << (32 - s));

    m_groups[m_groups.size() - 1] >>= s;

    emsmallen();

    return *this;
}

BigInt BigInt::operator<<(int rhs) const
{
    auto result = *this;

    if (rhs == 0)
        return result;

    if (rhs < 0)
        return result >>= -rhs;

    return result <<= rhs;
}

BigInt BigInt::operator+(BigInt const& rhs) const
{
    auto lhs = BigInt { *this };
    return lhs += rhs;
}

BigInt BigInt::operator-(BigInt const& rhs) const
{
    auto lhs = BigInt { *this };
    return lhs -= rhs;
}

BigInt BigInt::operator-() const
{
    auto rhs = BigInt { *this };
    rhs.m_negative ^= 1;
    return rhs;
}

bool BigInt::operator==(BigInt const& rhs) const
{
    if (m_negative != rhs.m_negative
        || m_groups.size() != rhs.m_groups.size()
        || size() != rhs.size())
        return false;

    for (auto lit = m_groups.begin(), rit = rhs.m_groups.begin(); lit != m_groups.end(); lit++, rit++)
        if (*lit != *rit)
            return false;

    return true;
}

bool BigInt::operator!=(BigInt const& rhs) const { return !(*this == rhs); }

int BigInt::operator<=>(BigInt const& rhs) const
{
    if (*this == rhs)
        return 0;

    if (m_negative && !rhs.m_negative)
        return 1;

    if (!m_negative && rhs.m_negative)
        return -1;

    if (m_groups.size() > rhs.m_groups.size() || size() > rhs.size())
        return (m_negative && rhs.m_negative) ? -1 : 1;

    if (m_groups.size() < rhs.m_groups.size() || size() < rhs.size())
        return (m_negative && rhs.m_negative) ? 1 : -1;

    for (auto lit = m_groups.begin(), rit = rhs.m_groups.begin(); lit != m_groups.end(); lit++, rit++) {
        if (*lit == *rit)
            continue;

        if (*lit < *rit)
            return (m_negative && rhs.m_negative) ? 1 : -1;
    }

    return 0;
}

bool BigInt::operator<=(BigInt const& rhs) const { return (*this <=> rhs) <= 0; }
bool BigInt::operator>=(BigInt const& rhs) const { return (*this <=> rhs) >= 0; }
bool BigInt::operator<(BigInt const& rhs) const { return (*this <=> rhs) < 0; }
bool BigInt::operator>(BigInt const& rhs) const { return (*this <=> rhs) > 0; }

std::ostream& operator<<(std::ostream& stream, BigInt const& number)
{
    auto const& group = number.m_groups;
    auto size = group.size();

    if (!size)
        return stream << 0;

    if (number.m_negative)
        stream << '-';

    if (size == 1)
        return stream << group[0];

    auto static muladd10 = [&](std::deque<uint32_t> x, uint64_t mul, uint64_t add) -> std::deque<uint32_t> {
        auto z = std::deque<uint32_t>(2 * number.m_groups.size() + 1);
        auto static constexpr base = BigInt::base;

        for (auto i = 0; i < x.size(); i++) {
            auto product = static_cast<uint64_t>(x[i]) * mul + z[i];

            z[i] = product % base;
            z[i + 1] = product / base;
        }

        auto carry = add;
        for (auto i = 0; i < z.size(); i++) {
            auto sum = static_cast<uint64_t>(z[i]) + carry;

            z[i] = sum % base;
            carry = sum / base;
        }

        return z;
    };

    auto z = std::deque<uint32_t> {};

    for (auto git = number.m_groups.rbegin(); git != number.m_groups.rend(); git++)
        z = muladd10(z, 1ull << 32, *git);

    while (z.back() == 0 && z.size() > 1)
        z.pop_back();

    auto i = z.size() - 1;

    while (z[i] == 0)
        if (__builtin_sub_overflow_p(i--, 1, i))
            return stream << 0;

    auto it = z.rbegin() + (z.size() - i - 1);

    stream << *it++;

    for (; it != z.rend(); it++)
        if (*it)
            stream << *it;

    return stream;
}
