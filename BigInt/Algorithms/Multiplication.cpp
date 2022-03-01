#include <BigInt/Algorithms/Algorithms.h>
#include <BigInt/BigInt.h>

#include <cstdint>
#include <deque>

template <typename T>
[[gnu::flatten]] BigInt naive_muladd(BigInt const& x, BigInt const& mul, BigInt const* add, T&& operation)
{
    /**
     * Impl of Algorithm M from The Art of Computer Programming
     * Vol. 2 Seminumerical Algorithms 3rd ed. pg. 268 by Donald Knuth.
     */

    auto result = std::deque<uint32_t>(x.groups() + mul.groups());
    auto carry = 0ull;

    // If add == 0 then Algorithm M degenerates into regular multiplication
    // It is also acceptible to pass nullptr instead of constructing a new zero valued BigInt
    if (add != nullptr)
        result.insert(result.begin(), add->get_groups().begin(), add->get_groups().end());

    for (auto i = 0uz; i < x.groups(); i++) {
        carry = 0;

        for (auto j = 0uz; j < mul.groups(); j++) {
            auto product = static_cast<uint64_t>(x.get_groups()[i]) * static_cast<uint64_t>(mul.get_groups()[j]) + result[i + j] + carry;

            // Specific operations for muladd dependent on base
            operation(product, carry, result[i + j]);
        }

        result[i + mul.groups()] = carry;
    }

    emsmallen(result);

    return { result };
}

[[gnu::flatten]] BigInt naive_muladd(BigInt const& x, BigInt const& mul, BigInt const* add)
{
    // naive_muladd for radix-b of 2^32

    return naive_muladd(x, mul, add,
                        [](auto const& product, auto& carry, auto& group) -> void {
                            group = static_cast<uint32_t>(product);
                            carry = product >> 32;
                        });
}

[[gnu::flatten]] BigInt multiply(BigInt const& x, BigInt const& y)
{
    // TODO: Implement Karatsuba and Toom-k algorithms for asymptotically faster results

    return naive_muladd(x, y, nullptr);
}
