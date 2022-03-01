#include <BigInt/Algorithms/Algorithms.h>
#include <BigInt/BigInt.h>

#include <stdexcept>

BigInt knuth(BigInt const& x, BigInt const& y, bool remainder)
{
    // Algorithm D, The Art of Computer Programming Vol. 2 Seminumerical Algorithms 3rd ed. pg. 272
    // Hacker's Delight divmnu64.c

    if (y.abs() == 0)
        throw new std::runtime_error("[BigInt] Div by 0.");

    if (y.groups() == 1)
        return knuth(x, y.get_groups()[0], remainder);

    if (x.abs() < y.abs()) {
        if (remainder)
            return x;

        return { 0 };
    }

    if (x.abs() == y.abs()) {
        if (remainder)
            return { 0 };

        return { 1 };
    }

    auto Q = std::deque<uint32_t>(x.groups());
    auto S = __builtin_clz(y.get_groups().back());
    auto U = (x << S).get_groups();
    auto V = (y << S).get_groups();

    emsmallen(U);
    emsmallen(V);

    U.push_back(0); // |U| = m + n + 1

    auto n = V.size();
    auto m = U.size() - n;

    for (auto j = m; j-- > 0;) {
        auto qhat = ((static_cast<uint64_t>(U[n + j]) << 32) | U[n + j - 1]) / V.back();
        auto rhat = ((static_cast<uint64_t>(U[n + j]) << 32) | U[n + j - 1]) % V.back();

        while (qhat >> 32 || qhat * V[n - 2] > ((rhat << 32) | U[n + j - 2])) {
            qhat--;
            rhat += V.back();

            if (rhat >> 32)
                break;
        }

        auto k = int64_t {};
        auto t = int64_t {};
        for (auto i = 0uz; i < n; i++) {
            auto p = qhat * V[i];
            t = static_cast<int64_t>(U[i + j]) - static_cast<uint32_t>(p) - k;
            U[i + j] = t;
            k = (p >> 32) - (t >> 32);
        }

        t = U[n + j] - k;
        U[n + j] = t;

        Q[j] = qhat;
        if (t < 0) {
            Q[j]--;
            k = 0;

            for (auto i = 0uz; i < n; i++) {
                t = static_cast<uint64_t>(U[i + j]) + V[i] + k;
                U[i + j] = t;
                k = t >> 32;
            }

            U[n + j] += k;
        }
    }

    if (remainder)
        return BigInt { U } >> S;

    emsmallen(Q);

    return { Q };
}

BigInt knuth(BigInt const& x, uint64_t y, bool remainder)
{
    // Degenerate case of Algorithm D
    auto Q = std::deque<uint32_t>(x.groups());
    auto k = uint64_t {};

    for (auto j = x.groups(); j-- > 0;) {
        uint64_t t = (k << 32) + x.get_groups()[j];
        Q[j] = t / y;
        k = t - Q[j] * y;
    }

    if (remainder)
        return BigInt(std::deque<uint32_t> { static_cast<uint32_t>(k), static_cast<uint32_t>(k >> 32) });

    return BigInt { Q };
}
