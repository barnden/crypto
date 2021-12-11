#pragma once

#include <cstdint>
#include <utility>

uint64_t gcd(uint64_t a, uint64_t b);

uint64_t Totient(uint64_t n);

std::pair<int64_t, int64_t> BezoutCoefficients(int64_t a, int64_t b);

uint64_t Modinv(uint64_t n, uint64_t mod);

uint64_t Modexp(uint64_t base, uint64_t exp, uint64_t mod);

bool MillerRabin(uint64_t n);

inline uint64_t Modsub(uint64_t a, uint64_t b, uint64_t mod)
{
    // Compute (a - b (mod m)) (mod 2^64)
    if (a == b)
        return 0;

    while (a < b)
        a += mod;

    return (a - b) % mod;
}
