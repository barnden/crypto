#include <BigInt/BigInt.h>
#include <EllipticCurve/EllipticCurve.h>
#include <Modmath.h>

#include <iostream>
#include <random>

BigInt gcd(BigInt const& a, BigInt const& b)
{
    // TODO: Implement binary gcd for faster (?) computations
    // On my system I've found binary gcd to be slower than just Euclid's
    // for relatively small (~2bn) integers. Maybe something wrong with my impl?
    if (a == 0)
        return b.abs();

    return gcd(b.abs() % a.abs(), a.abs());
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

std::pair<BigInt, BigInt> BezoutCoefficients(BigInt a, BigInt b)
{
    // Compute Bezout Coeffs. s, t for a, b such that s * a + t * b = (a, b)
    if (b < a) {
        auto c = b;
        b = a;
        a = c;
    }

    auto pr = a;
    auto r = b;

    auto ps = BigInt { 1 };
    auto s = BigInt {};

    auto pt = BigInt {};
    auto t = BigInt { 1 };

    while (r != 0) {
        auto q = pr / r;
        auto temp = BigInt {};

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

BigInt Modinv(BigInt const& n, BigInt const& mod)
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
    if (exp > mod)
        exp %= Totient(mod);

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

    if (n.get_groups().front() % 2 == 0)
        return true;

    auto np = n - 1;
    auto r = np.trailing_zeros();
    auto d = np >> r;

    // clang-format off
    auto static constexpr bases = std::integer_sequence<uint64_t, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43> {};
    auto static constexpr large_bases = std::integer_sequence<uint64_t, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
    367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
    503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647,
    653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971,
    977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097,
    1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237,
    1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399,
    1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523,
    1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
    1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801,
    1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973,
    1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029> {};
    // clang-format on

    bool composite = false;

    auto static MRTest = [&](BigInt const& base) {
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

    [&]<std::size_t... I>(std::integer_sequence<uint64_t, I...>) { (MRTest(I), ...); }
    (large_bases);

    // Otherwise, we perform MR using randomized bases
    BigInt base {};
    for (auto i = 0; i < 10; i++) {
        if (composite)
            return true;

        base.random(n.size() - 1);

        MRTest(base);
    }

    return composite;
}

BigInt LenstraFactorization(BigInt const& n)
{
    // Find a nontrival factor of n using Lenstra's factorization method
    // Works best for n semiprime, i.e. n = pq where p and q distinct primes and q of much smaller order than p

    // FIXME: Make this interface nicer, perhaps static random to construct new BigInt directly
    // Something like: auto a  = BigInt.random() % n;

    auto a = BigInt();
    auto x = BigInt();
    auto y = BigInt();

    a.random(n.size());
    x.random(n.size());
    y.random(n.size());

    a %= n;
    x %= n;
    y %= n;

    // FIXME: a - b - c != a - (b + c)
    // Currently: a - (b + c) gives the expected result for a - b - c, which is what we use below.
    auto b = ((y * y) - ((x * x * x) + (x * a))) % n;

    auto ec = EllipticCurve(a, b, n);
    auto P = Point(x, y, ec);

    auto j = 1;

    while (j++) {
        auto Q = j * P;

        if (!Q.get_w()) {
            auto d = gcd(Q.get_x() - P.get_x(), n);

            if (d == 1 || d == n)
                return LenstraFactorization(n);

            return d;
        }

        P = Q;
    }

    // This is here to make the compiler happy
    return 0;
}
