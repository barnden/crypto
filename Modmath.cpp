#include "Modmath.h"
#include <random>

BigInt gcd(BigInt const& a, BigInt const& b)
{
    // TODO: Implement binary gcd for faster (?) computations
    // On my system I've found binary gcd to be slower than just Euclid's
    // for relatively small (~2bn) integers. Maybe something wrong with my impl?
    if (a == 0)
        return b;

    return gcd(b % a, a);
}

BigInt Totient(BigInt const& n)
{
    // Return Euler Totient Function of n
    auto accumulator = BigInt { 1 };

    if (n == 0)
        return { 0 };

    if (n == 1)
        return accumulator;

    // Totient(p) = p - 1, where p prime.
    // Check if number is prime
    if (n % 2 == 1 && !MillerRabin(n))
        return n - 1;

    /** Optimization:
     * Totient(p^q) = (p - 1)p^(q - 1)
     * For p = 2, Totient(2^q) = 2^(q-1)
     * Using multiplicity of Euler's Totient Function:
     * Totient(2^q * r)
     *  = Totient(2^q) * Totient(r)
     *  = 2^(q-1) * Totient(r)
     *  = (1 << q - 1) * Totient(r)
     */

    auto k = n.trailing_zeros();

    if (k)
        return Totient(n >> k) << (k - 1);

    // FIXME: This part is exponential, we can achieve sub-exponential results using Lenstra/Quadratic Sieve/GNFS methods.
    for (BigInt i = 3; i < n; i += 1)
        if (gcd(i, n) == 1)
            accumulator += 1;

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

BigInt Modexp(BigInt const& base, BigInt exp, BigInt const& mod)
{
    // Implementation of fast powering approach to exponentiation in integer rings
    // TODO: Handle negative exponents, i.e. find d := Modinv(base, mod), then return Modexp(d, exp, mod).

    // If a = 0 (mod m), then a^n = 0 (mod m) for all n in Z/mZ
    if (base % mod == 0)
        return 0;

    // Short circuit squaring
    if (exp == 2)
        return (base * base) % mod;

    // Given a^n (mod m); If n > m, then exp := exp (mod Totient(m))
    // if (exp > mod)
    //     exp %= Totient(mod);

    auto accumulator = BigInt { 1 };

#ifdef MONTGOMERY
    // Use Montgomery ladder method to do fast powering
    auto g = base;
    for (auto i = exp.size(); i-- > 0;) {
        if (exp.bit_at(i)) {
            accumulator = (accumulator * g) % mod;
            g = (g * g) % mod;
            continue;
        }

        g = (accumulator * g) % mod;
        accumulator = (accumulator * accumulator) % mod;
    }
#else
    // Use traditional fast powering; suceptible to side channel attacks
    for (auto i = exp.size(); i-- > 0;) {
        accumulator = (accumulator * accumulator) % mod;

        if (exp.bit_at(i))
            accumulator = (accumulator * base) % mod;
    }
#endif

    return accumulator % mod;
}

bool MillerRabin(BigInt const& n)
{
    /**
     * Miller-Rabin primality test
     * Returns true if composite; false if probably prime
     *
     * Miller-Rabin is a non-deterministic primality test.
     * It has been shown by Sorenson and Webster (doi:10.1090/mcom/3134) that for
     * a composite number n < 3,317,044,064,679,887,385,961,981 (< 82 bits) at
     * least one of the bases below will be a witness to the compositeness of n.
     */

    if (n.get_groups()[0] % 2 == 0)
        return true;

    auto np = n - 1;
    auto r = np.trailing_zeros();
    auto d = np >> r;

    auto static constexpr bases = std::integer_sequence<uint64_t, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43> {};
    bool composite = false;

    auto static MRTest = [&](std::size_t base) {
        if (base == n)
            return;

        if (composite)
            return;

        auto x = Modexp(base, d, n);

        if (x == 1 || x == np)
            return;

        for (auto i = 0ull; i < r; i++) {
            x = (x * x) % n;

            if (x == np)
                return;
        }

        composite = true;
    };

    [&]<std::size_t... I>(std::integer_sequence<uint64_t, I...>) { (MRTest(I), ...); }
    (bases);

    // Early return if bit length is less than 82; see comments above.
    if (n.size() < 82 || composite)
        return composite;

    // Otherwise, we perform MR using randomized bases
    auto static rd = std::random_device {};
    auto static e2 = std::mt19937_64 { rd() };
    auto static dist = std::uniform_int_distribution<uint64_t> { 1ull << 61, 1ull << 62 };

    for (auto i = 0; i < 40; i++) {
        if (composite)
            return true;

        MRTest(dist(e2));
    }

    return composite;
}
