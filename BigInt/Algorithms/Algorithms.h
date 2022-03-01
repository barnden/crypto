#pragma once

#include <cstdint>

class BigInt;

// Multiplication Algorithms
template <typename T>
BigInt naive_muladd(BigInt const& x, BigInt const& mul, BigInt const* add, T&& operation);
BigInt naive_muladd(BigInt const& x, BigInt const& mul, BigInt const* add);
BigInt multiply(BigInt const& x, BigInt const& y);

// Division Algorithms
BigInt knuth(BigInt const& x, uint64_t y, bool remainder);
BigInt knuth(BigInt const& x, BigInt const& y, bool remainder);
