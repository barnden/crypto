#pragma once

#include "BigInt.h"
#include <cstdint>
#include <utility>

BigInt gcd(BigInt const& a, BigInt const& b);

BigInt Totient(BigInt const& n);

std::pair<int64_t, int64_t> BezoutCoefficients(int64_t a, int64_t b);

uint64_t Modinv(uint64_t n, uint64_t mod);

BigInt Modexp(BigInt const& base, BigInt exp, BigInt const& mod);

bool MillerRabin(BigInt const& n);

inline uint64_t Modsub(uint64_t a, uint64_t b, uint64_t mod)
{
    // Compute (a - b (mod m)) (mod 2^64)
    if (a == b)
        return 0;

    while (a < b)
        a += mod;

    return (a - b) % mod;
}
