#pragma once

#include <cstdint>
#include <deque>
#include <string>
#include <utility>

#include <BigInt/Algorithms/Algorithms.h>

void emsmallen(std::deque<uint32_t>& groups);

template <typename T>
concept Numeric = std::convertible_to<T, std::size_t>;

template <Numeric T, size_t R>
size_t constexpr get_max_digits()
{
    if (R < 2)
        return 0;

    size_t digits = 0;

    T radix = R;
    T place = 1;

    do {
        place *= radix;
        digits++;
    } while (!__builtin_mul_overflow_p(place, radix, place));

    return digits;
}

template <Numeric T, size_t R>
T constexpr get_base()
{
    auto constexpr max_digits = get_max_digits<T, R>();

    T base = 1;
    for (auto i = 0uz; i < max_digits; i++)
        base *= R;

    return base;
}

class BigInt {

private:
    std::deque<uint32_t> m_groups;
    bool m_negative;

    void embiggen(BigInt const& other);
    void embiggen(size_t size);
    void emsmallen();

    // TODO: Make this work for radices not 10
    size_t static constexpr radix = 10;
    size_t static constexpr digits = get_max_digits<uint32_t, radix>();
    size_t static constexpr base = get_base<uint32_t, radix>();
    size_t static constexpr base_sz = sizeof(uint32_t) * 8;

public:
    BigInt();
    BigInt(int64_t value);
    BigInt(std::string number);
    BigInt(std::deque<uint32_t> group);

    BigInt& operator+=(BigInt const& rhs);
    BigInt& operator-=(BigInt const& rhs);
    BigInt& operator*=(BigInt const& lhs);
    BigInt& operator/=(BigInt const& rhs);
    BigInt& operator%=(BigInt const& rhs);

    BigInt& operator*=(int rhs);
    BigInt& operator/=(int rhs);
    BigInt& operator%=(int rhs);

    BigInt& operator*=(uint64_t rhs);
    BigInt& operator/=(uint64_t rhs);
    BigInt& operator%=(uint64_t rhs);

    BigInt operator+(BigInt const& rhs) const;
    BigInt operator-(BigInt const& rhs) const;
    BigInt operator-() const;
    BigInt operator*(BigInt const& rhs) const;
    BigInt operator/(BigInt const& rhs) const;
    BigInt operator%(BigInt const& rhs) const;

    BigInt operator*(int rhs) const;
    BigInt operator/(int rhs) const;
    BigInt operator%(int rhs) const;

    BigInt operator*(uint64_t rhs) const;
    BigInt operator/(uint64_t rhs) const;
    BigInt operator%(uint64_t rhs) const;

    bool operator==(BigInt const& rhs) const;
    bool operator!=(BigInt const& rhs) const;
    int operator<=>(BigInt const& rhs) const;
    bool operator<=(BigInt const& rhs) const;
    bool operator>=(BigInt const& rhs) const;
    bool operator<(BigInt const& rhs) const;
    bool operator>(BigInt const& rhs) const;

    BigInt& operator<<=(int rhs);
    BigInt& operator>>=(int rhs);

    BigInt operator<<(int rhs) const;
    BigInt operator>>(int rhs) const;

    size_t size() const;
    inline size_t groups() const { return m_groups.size(); };

    friend std::ostream& operator<<(std::ostream& stream, BigInt const& number);

    inline std::deque<uint32_t> const& get_groups() const { return m_groups; }
    inline bool is_negative() const { return m_negative; }
    size_t trailing_zeros() const;
    bool bit_at(size_t n) const;
    BigInt abs() const;
    void random(int bits);
    bool is_power_of_two() const;
};
