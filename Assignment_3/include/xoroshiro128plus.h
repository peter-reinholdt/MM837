#ifndef XORO_HEADER
#define XORO_HEADER

#include <stdint.h>
#include <array>
#include <random>

struct xoroshiro128plus{
    private:
        uint64_t state[2];

        static inline uint64_t rotl(const uint64_t x, int k) {
            return (x << k) | (x >> (64 - k));
        }
    public:
        using result_type = uint64_t;
        constexpr static result_type min() {return 0;}
        constexpr static result_type max() {return -1;}
        result_type operator()();
        xoroshiro128plus(const std::array<uint32_t, 4> &);
};
#endif
