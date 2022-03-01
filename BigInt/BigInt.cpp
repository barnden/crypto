#include <BigInt/Algorithms/Algorithms.h>
#include <BigInt/BigInt.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>

BigInt::BigInt()
    : m_groups({})
    , m_negative(false)
{
    m_groups.push_back(0);
}

BigInt::BigInt(int64_t number)
    : m_groups({})
    , m_negative(false)
{
    if (number < 0ll) {
        m_negative = true;
        number *= -1;
    }

    auto static constexpr offset = sizeof(uint32_t) * 8;

    if (number == 0) {
        m_groups.push_back(0);
    } else {
        while (number) {
            m_groups.push_back(static_cast<uint32_t>(number));

            number >>= offset;
        }
    }
}

BigInt::BigInt(std::string number)
    : m_groups({})
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
    auto add = BigInt {};

    for (auto git = m_groups.rbegin(); git != m_groups.rend(); git++) {
        add.m_groups[0] = *git;
        z = naive_muladd(z, base, &add);
    }

    m_groups = z.m_groups;

    emsmallen();
}

BigInt::BigInt(std::deque<uint32_t> group)
    : m_groups(group)
    , m_negative(false)
{
    emsmallen();
}

void BigInt::random(int bits)
{
    // Generate a random big integer
    // Utilizes a Mersenne twister -> determinstic, not suitable for actual cryptography
    auto static rd = std::random_device {};
    auto static e2 = std::mt19937_64 { rd() };
    auto static dist = std::uniform_int_distribution<uint64_t> { 0, 1ull << 32 };
    m_groups.clear();

    for (auto i = 0; i < bits / 32; i++)
        m_groups.push_back(static_cast<uint32_t>(dist(e2)));

    m_groups.push_back(static_cast<uint32_t>(dist(e2)) >> (32 - (bits % 32)));
}

size_t BigInt::trailing_zeros() const
{
    auto i = 0uz;

    for (; i < m_groups.size(); i++)
        if (m_groups[i] != 0)
            break;

    return (32 * i) + __builtin_ctz(m_groups[i]);
}

size_t BigInt::size() const
{
    // Get size of number in bits
    return ((groups() - 1) * 32) + (32 - __builtin_clz(m_groups.back()));
}

bool BigInt::is_power_of_two() const
{
    return trailing_zeros() == size() - 1;
}

bool BigInt::bit_at(size_t n) const
{
    auto const& group = m_groups[n / 32];

    return group & (1ul << (n % 32));
}

BigInt BigInt::abs() const
{
    auto ret = *this;

    ret.m_negative = false;

    return ret;
}

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

    for (auto i = 0uz; i <= size - groups(); i++)
        m_groups.push_back(0);
}

inline void emsmallen(std::deque<uint32_t>& groups)
{
    // Remove all leading digit groups valued at 0, except for the last one
    if (groups.size() == 0) {
        // This case should never happen
        groups.push_back(0);
        return;
    }

    while (groups.back() == 0 && groups.size() > 1)
        groups.pop_back();
}

inline void BigInt::emsmallen() { ::emsmallen(m_groups); }

BigInt& BigInt::operator-=(BigInt const& rhs)
{
    if (rhs.m_negative)
        return *this += -rhs;

    auto it = m_groups.begin();

    auto lit = std::as_const(m_groups).begin();
    auto lend = std::as_const(m_groups).end();

    auto rit = rhs.m_groups.begin();
    auto rend = rhs.m_groups.end();

    if (*this < rhs) {
        std::swap(lit, rit);
        std::swap(lend, rend);

        embiggen(rhs);
        rend = std::as_const(m_groups).end();

        m_negative = true;
    }

    auto borrow = int64_t {};

    for (; lit != lend; lit++, it++) {
        auto difference = static_cast<int64_t>(*lit) + borrow;

        if (rit != rend) {
            difference -= static_cast<int64_t>(*rit);
            rit++;
        }

        *it = static_cast<uint32_t>(difference);
        borrow = difference >> 32;
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

    m_negative ^= rhs.m_negative;
    m_groups = multiply(*this, rhs).m_groups;

    emsmallen();

    return *this;
}

BigInt& BigInt::operator*=(int rhs)
{
    if (rhs < 0) {
        m_negative ^= 1;
        rhs *= -1;
    }

    return *this *= static_cast<uint64_t>(rhs);
}

BigInt& BigInt::operator*=(uint64_t rhs)
{
    m_groups = multiply(*this, rhs).m_groups;

    emsmallen();

    return *this;
}

BigInt BigInt::operator*(BigInt const& rhs) const { return BigInt { *this } *= rhs; }
BigInt BigInt::operator*(int rhs) const { return BigInt { *this } *= rhs; }
BigInt BigInt::operator*(uint64_t rhs) const { return BigInt { *this } *= rhs; }

BigInt& BigInt::operator/=(BigInt const& rhs)
{
    m_negative ^= rhs.m_negative;
    m_groups = knuth(*this, rhs, false).m_groups;

    emsmallen();

    return *this;
}

BigInt& BigInt::operator/=(int rhs)
{
    if (rhs < 0) {
        m_negative ^= 1;
        rhs *= -1;
    }

    return *this /= static_cast<uint64_t>(rhs);
}

BigInt& BigInt::operator/=(uint64_t rhs)
{
    m_groups = knuth(*this, rhs, false).m_groups;

    emsmallen();

    return *this;
}

BigInt BigInt::operator/(BigInt const& rhs) const { return BigInt { *this } /= rhs; }
BigInt BigInt::operator/(int rhs) const { return BigInt { *this } /= rhs; }
BigInt BigInt::operator/(uint64_t rhs) const { return BigInt { *this } /= rhs; }

BigInt& BigInt::operator%=(BigInt const& rhs)
{
    if (rhs.m_negative)
        throw new std::runtime_error("[BigInt] Negative modulus");

    m_groups = knuth(*this, rhs, true).m_groups;

    while (m_negative)
        *this += rhs;

    emsmallen();

    return *this;
}

BigInt& BigInt::operator%=(int rhs)
{
    if (rhs < 0)
        throw new std::runtime_error("[BigInt] Negative modulus");

    return *this %= static_cast<uint64_t>(rhs);
}

BigInt& BigInt::operator%=(uint64_t rhs)
{
    m_groups = knuth(*this, rhs, true).m_groups;

    if (m_negative)
        *this += rhs;

    emsmallen();

    return *this;
}

BigInt BigInt::operator%(BigInt const& rhs) const { return BigInt { *this } %= rhs; }
BigInt BigInt::operator%(int rhs) const { return BigInt { *this } %= rhs; }
BigInt BigInt::operator%(uint64_t rhs) const { return BigInt { *this } %= rhs; }

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

    if (static_cast<size_t>(rhs) >= size()) {
        m_groups.clear();
        m_groups[0] = 0;

        return *this;
    }

    auto groups = static_cast<size_t>(rhs / 32);

    if (groups > m_groups.size()) {
        m_groups.clear();
        m_groups.push_back(0);

        return *this;
    }

    for (auto i = 0uz; i < groups; i++)
        m_groups.pop_front();

    auto s = rhs % 32;
    for (auto i = 0uz; i < m_groups.size() - 1; i++)
        m_groups[i] = (m_groups[i] >> s) | static_cast<uint64_t>(m_groups[i + 1] << (32 - s));

    m_groups[m_groups.size() - 1] >>= s;

    emsmallen();

    return *this;
}

BigInt BigInt::operator>>(int rhs) const
{
    auto result = *this;

    if (rhs == 0)
        return result;

    if (rhs < 0)
        return result <<= -rhs;

    return result >>= rhs;
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
        auto z = std::deque<uint32_t>(4 * (number.groups() + 1));
        auto static constexpr base = BigInt::base;

        for (auto i = 0uz; i < x.size(); i++) {
            auto product = static_cast<uint64_t>(x[i]) * mul + z[i];

            z[i] = product % base;
            z[i + 1] = product / base;
        }

        auto carry = add;
        for (auto i = 0uz; i < z.size(); i++) {
            auto sum = static_cast<uint64_t>(z[i]) + carry;

            z[i] = sum % base;
            carry = sum / base;
        }

        return z;
    };

    auto z = std::deque<uint32_t> {};

    for (auto git = number.m_groups.rbegin(); git != number.m_groups.rend(); git++)
        z = muladd10(z, 1ull << 32, *git);

    auto i = z.size() - 1;

    while (z[i] == 0)
        if (__builtin_sub_overflow_p(i - 1, 1, i))
            return stream << 0;

    auto it = z.rbegin() + (z.size() - i - 1);

    stream << *it++;

    for (; it != z.rend(); it++)
        if (*it)
            stream << *it;

    return stream;
}
