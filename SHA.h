#pragma once
#include "Modmath.h"
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace SHA256 {
inline uint32_t ROTR(uint32_t x, int k)
{
    return (x >> k) | (x << (32 - k));
}

std::array<uint32_t, 8> Hash(std::vector<bool> message)
{
    uint32_t static constexpr k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };

    auto static h = std::array<uint32_t, 8> { 0x6a09e667,
                                              0xbb67ae85,
                                              0x3c6ef372,
                                              0xa54ff53a,
                                              0x510e527f,
                                              0x9b05688c,
                                              0x1f83d9ab,
                                              0x5be0cd19 };

    auto static a = h;

    auto const L = message.size();
    auto const K = Modsub(512, L + 65, 512);
    auto padding = std::vector<bool>(K + 65);

    for (auto i = 0; i < 64; i++)
        padding[K + 1 + i] = (L >> (63 - i)) & 1;

    message.reserve(L + K + 65);
    message.insert(message.end(), padding.begin(), padding.end());

    message[L] = true;

    assert(message.size() % 512 == 0);

    for (auto i = 0; i < (message.size() >> 9); i++) {
        uint32_t w[64] = { 0 };

        // Copy first 512 bits of message into w[0..16]
        for (auto j = 0; j < 16; j++)
            for (auto k = 0; k < 32; k++)
                w[j] |= message[j * 32 + k] << (31 - k);

        for (auto j = 16; j < 64; j++)
            w[j] = w[j - 16] + w[j - 7]
                   + (ROTR(w[j - 15], 7) ^ ROTR(w[j - 15], 18) ^ (w[j - 15] >> 3))
                   + (ROTR(w[j - 2], 17) ^ ROTR(w[j - 2], 19) ^ (w[j - 2] >> 10));

        for (auto j = 0; j < 64; j++) {
            auto temp1 = a[7] + k[j] + w[j]
                         + (a[6] ^ (a[4] & (a[5] ^ a[6])))
                         + (ROTR(a[4], 6) ^ ROTR(a[4], 11) ^ ROTR(a[4], 25));
            auto temp2 = ((ROTR(a[0], 2) ^ ROTR(a[0], 13) ^ ROTR(a[0], 22)))
                         + ((a[0] & a[1]) | (a[2] & (a[0] | a[1])));

            for (auto i = 7; i-- > 0;)
                a[i + 1] = a[i];

            a[4] += temp1;
            a[0] = temp1 + temp2;
        }

        [&]<std::size_t... I>(std::index_sequence<I...>)
        {
            ((h[I] += a[I]), ...);
        }
        (std::make_index_sequence<8> {});
    }

    return h;
}
}
