#include "BigInt.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

BigInt::BigInt()
    : m_negative(false)
{
}

BigInt::BigInt(uint64_t number)
    : m_negative(false)
{
    if (number < 0) {
        m_negative = true;
        number *= -1;
    }

    while (number) {
        m_group.push_back(number % base);

        number /= base;
    }
}

BigInt::BigInt(std::string number)
{
    m_group.reserve((number.size() / digits) + 1);
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
        auto num = base_t {};
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

        if (num || !m_group.size())
            m_group.push_back(num);
    }
}

void BigInt::embiggen(BigInt const& other)
{
    // Pad m_group to be other.m_group.size + 1
    if (m_group.size() > other.m_group.size())
        return;

    embiggen(other.m_group.size() + 1);
}

void BigInt::embiggen(size_t size)
{
    if (m_group.size() > size)
        return;

    auto sz = m_group.size();

    m_group.reserve(size);

    for (auto i = 0; i <= size - sz; i++)
        m_group.push_back(0);
}

void BigInt::emsmallen()
{
    // Remove all leading digit groups valued at 0, except for the last one
    while (m_group.back() == 0 && m_group.size() > 1)
        m_group.pop_back();
}

size_t BigInt::size() const
{
    size_t n = digits * (m_group.size() - 1);

    for (auto i = m_group.back(); i > 0; i /= 10)
        n++;

    return n;
}

BigInt& BigInt::operator-=(BigInt const& rhs)
{
    if (rhs.m_negative)
        return *this += -rhs;

    embiggen(rhs);

    auto lit = m_group.begin();
    auto const& lend = m_group.end();

    auto rit = std::as_const(rhs.m_group).begin();
    auto const& rend = rhs.m_group.end();

    auto sum = acc_t {};
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
            *std::prev(lit) = sum % base;
            sum /= base;

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

    auto lit = m_group.begin();
    auto const& lend = m_group.end();

    auto rit = std::as_const(rhs.m_group).begin();
    auto const& rend = rhs.m_group.end();

    auto sum = acc_t {};
    while (lit != lend || rit != rend) {
        if (lit != lend) {
            sum += *lit;
            lit++;
        }

        if (rit != rend) {
            sum += *rit;
            rit++;
        }

        *std::prev(lit) = sum % base;
        sum /= base;
    }

    if (sum)
        m_group.push_back(1);

    emsmallen();

    return *this;
}

void naive_multiplication(BigInt& a, BigInt const& b)
{
    using acc_t = BigInt::acc_t;
    using base_t = BigInt::base_t;
    size_t static const constexpr offset = sizeof(base_t) * 8;

    auto const& x = a.m_group;
    auto const& y = b.m_group;

    auto z = std::vector<base_t>(x.size() + y.size() + 1);
    auto carry = acc_t {};

    for (auto i = 0; i < x.size(); i++) {
        carry = 0;

        for (auto j = 0, k = i; j < y.size(); j++, k++) {
            auto product = static_cast<acc_t>(x[i]) * static_cast<acc_t>(y[j]) + z[k] + carry;

            z[k] = product % BigInt::base;
            carry = product / BigInt::base;
        }

        z[i + y.size()] = carry;
    }

    a.m_group = z;
    a.emsmallen();
}

BigInt& BigInt::operator*=(BigInt const& rhs)
{
    // Perform multiplication
    // TODO: Implement Karatsuba and Toom-k

    m_negative ^= rhs.m_negative;
    naive_multiplication(*this, rhs);

    return *this;
}

BigInt& BigInt::operator<<=(int rhs)
{
    if (rhs == 0)
        return *this;

    if (rhs < 0)
        return *this >>= -rhs;

    // TODO: Implement then use multiplication by integer primitive
    BigInt r;

    if (rhs < 63)
        r = { 1ull << rhs };
    else {
        // TODO: There has to be a smarter way
        r = { 1ull };

        for (auto i = 0; i < rhs; i++)
            r += r;
    }

    return *this *= r;
}

BigInt& BigInt::operator>>=(int rhs)
{
    // TODO: Study Knuth's division algorithm

    if (rhs == 0)
        return *this;

    if (rhs < 0)
        return *this <<= -rhs;

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
        || m_group.size() != rhs.m_group.size()
        || size() != rhs.size())
        return false;

    for (auto lit = m_group.begin(), rit = rhs.m_group.begin(); lit != m_group.end(); lit++, rit++)
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

    if (m_group.size() > rhs.m_group.size() || size() > rhs.size())
        return (m_negative && rhs.m_negative) ? -1 : 1;

    if (m_group.size() < rhs.m_group.size() || size() < rhs.size())
        return (m_negative && rhs.m_negative) ? 1 : -1;

    for (auto lit = m_group.begin(), rit = rhs.m_group.begin(); lit != m_group.end(); lit++, rit++) {
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
    auto const& group = number.m_group;
    auto size = group.size();

    if (!size)
        return stream << 0;

    if (number.m_negative)
        stream << '-';

    if (size == 1)
        return stream << group[0];

    auto i = size - 1;

    while (group[i] == 0)
        if (__builtin_sub_overflow_p(i--, 1, i))
            return stream << 0;

    auto it = group.rbegin() + (size - i - 1);

    stream << *it++;

    for (; it != group.rend(); it++)
        if (*it)
            stream << *it;

    return stream;
}
