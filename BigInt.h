#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

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
    for (auto i = 0; i < max_digits; i++)
        base *= R;

    return base;
}

class BigInt {
private:
    using base_t = uint32_t;
    using acc_t = uint64_t;

    std::vector<base_t> m_group;
    bool m_negative;

    void embiggen(BigInt const& other);
    void embiggen(size_t size);
    void emsmallen();

    // TODO: Make this work for radices not 10
    size_t static constexpr radix = 10;
    size_t static constexpr digits = get_max_digits<base_t, radix>();
    size_t static constexpr base = get_base<base_t, radix>();
    size_t static constexpr base_sz = sizeof(base_t) * 8;

    // Multiplication algorithms
    template <Numeric T, Numeric U> friend std::vector<T> naive_multiplication(std::vector<T> const& x, U y); // O(n^2)
    template <Numeric T, Numeric U> friend std::vector<T> naive_multiplication(std::vector<T> const& x, std::vector<T> const& y); // O(n^2)
    template <Numeric T, Numeric U> friend std::vector<T> naive_muladd(std::vector<T> const& x, U mul, U add); // O(n^2)
    friend BigInt karatsuba(BigInt const& a, BigInt const& b); // O(n^1.58)

    // Division algorithms
    template <Numeric T, Numeric U> friend std::vector<T> knuth(std::vector<T> const& x, U y); // O(n^2)
    template <Numeric T, Numeric U> friend std::vector<T> knuth(std::vector<T> const& x, std::vector<T> const& y); // O(n^2)
    friend BigInt knuth(BigInt const& a, BigInt const& b); // O(n^2)

public:
    BigInt();
    BigInt(uint64_t value);
    BigInt(std::string number);
    BigInt(std::vector<base_t> group);

    BigInt& operator+=(BigInt const& rhs);
    BigInt& operator-=(BigInt const& rhs);
    BigInt& operator*=(BigInt const& lhs);
    BigInt& operator/=(BigInt const& rhs);

    BigInt& operator*=(int rhs);

    BigInt operator+(BigInt const& rhs) const;
    BigInt operator-(BigInt const& rhs) const;
    BigInt operator-() const;

    bool operator==(BigInt const& rhs) const;
    bool operator!=(BigInt const& rhs) const;
    int operator<=>(BigInt const& rhs) const;
    bool operator<=(BigInt const& rhs) const;
    bool operator>=(BigInt const& rhs) const;
    bool operator<(BigInt const& rhs) const;
    bool operator>(BigInt const& rhs) const;

    // TODO: Implement bitshift operators
    BigInt& operator<<=(int rhs);
    BigInt& operator>>=(int rhs);

    BigInt operator<<(int rhs) const;
    BigInt operator>>(int rhs) const;

    size_t size() const;

    friend std::ostream& operator<<(std::ostream& stream, BigInt const& number);

    inline std::vector<base_t> const& groups() const { return m_group; }
    inline bool is_negative() const { return m_negative; }
};
