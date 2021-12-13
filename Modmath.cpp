#include "Modmath.h"
#include <algorithm>

uint64_t gcd(uint64_t a, uint64_t b)
{
    // Return (a, b); assume b > a
    if (a == 0)
        return b;

    return gcd(b % a, a);
}

uint64_t Totient(uint64_t n)
{
    // Return Euler Totient Function of n
    auto accumulator = 1ull;

    // Totient(p) = p - 1, where p prime.
    // Check if number is prime
    if (n > 41 && n % 2 && !MillerRabin(n))
        return n - 1;

    for (auto i = 2; i < n; i++)
        if (gcd(i, n) == 1)
            accumulator++;

    return accumulator;
}

std::pair<int64_t, int64_t> BezoutCoefficients(int64_t a, int64_t b)
{
    // Compute Bezout Coeffs. s, t for a, b such that s * a + t * b = (a, b)
    if (b < a) {
        auto c = b;
        b = a;
        a = c;
    }

    auto pr = a;
    auto r = b;

    auto ps = 1;
    auto s = 0;

    auto pt = 0;
    auto t = 1;

    while (r) {
        auto q = pr / r;
        auto temp = 0;

        temp = r;
        r = pr - q * r;
        pr = temp;

        temp = s;
        s = ps - q * s;
        ps = temp;

        temp = t;
        t = pt - q * t;
        pt = temp;
    }

    return std::make_pair(ps, pt);
}

uint64_t Modinv(uint64_t n, uint64_t mod)
{
    auto coeff = BezoutCoefficients(n % mod, mod).first;

    if (coeff < 0)
        return (coeff + mod) % mod;

    return coeff;
}

uint64_t Modexp(uint64_t base, uint64_t exp, uint64_t mod)
{
    // Implementation of fast powering approach to exponentiation in integer rings

    /**
     * Allow negative exponents?
     * Currently we assume positive exponent.
     * However, we could change exp to int64_t, then compute
     * a = Modexp(base, |exp|, mod), then return Modinv(a, mod)
     * if exp < 0.
     *
     * TODO: Verify this assumption:
     * Since uint64_t is essentially the same as (mod 2^64) (**).
     * Switching exp to int64_t won't change too much, i.e.
     * Modexp(2, 2^64 - 1, 97) in is the same as Modexp(2, -1, 97).
     *
     * (**): Technically overflow is undefined behavior in ISO C++.
     */

    // If a = 0 (mod m), then a^n = 0 (mod m) for all n in Z/mZ
    if (base % mod == 0)
        return 0;

    // Short circuit squaring
    if (exp == 2)
        return (base * base) % mod;

    // Given a^n (mod m); If n > m, then exp := exp (mod Totient(m))
    if (exp > mod)
        exp %= Totient(mod);

    auto const m = (sizeof(exp) * 8) - __builtin_clzll(exp);
    auto accumulator = 1ull;

#ifdef MONTGOMERY
    // Use Montgomery ladder method to do fast powering
    auto g = base;
    for (auto i = m; i-- > 0;) {
        if (exp & (1 << i)) {
            accumulator = (accumulator * g) % mod;
            g = (g * g) % mod;
            continue;
        }

        g = (accumulator * g) % mod;
        accumulator = (accumulator * accumulator) % mod;
    }
#else
    // Use traditional fast powering; suceptible to side channel attacks
    for (auto i = m; i-- > 0;) {
        accumulator = (accumulator * accumulator) % mod;

        if (exp & (1 << i))
            accumulator = (accumulator * base) % mod;
    }
#endif

    return accumulator;
}

bool MillerRabin(uint64_t n)
{
    /**
     * Miller-Rabin primality test
     * Returns true if composite; false if probably prime
     *
     * Miller-Rabin is a non-deterministic primality test.
     * It has been shown by Sorenson and Webster (doi:10.1090/mcom/3134) that for
     * a composite number n < 3,317,044,064,679,887,385,961,981 at least one of
     * the bases below will be a witness to the compositeness of n. Thus, for the
     * values of n less than the bound, Miller-Rabin is effectively deterministic.
     */

    if (n % 2 == 0)
        return true;

    auto r = 0;
    auto np = n - 1;
    auto d = np;

    while (d % 2 == 0) {
        d >>= 1;
        r++;
    }

    auto static constexpr bases = std::integer_sequence<uint64_t, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41> {};
    bool composite = false;

    [&]<std::size_t... I>(std::index_sequence<I...>)
    {
        ([&](std::size_t base) {
            if (composite)
                return;

            auto x = Modexp(base, d, n);

            if (x == 1 || x == np)
                return;

            for (auto i = 0ull; i < r; i++) {
                x = Modexp(base, 2, n);

                if (x == np)
                    return;
            }

            composite = true;
        }(I),
         ...);
    }
    (bases);

    return composite;
}
