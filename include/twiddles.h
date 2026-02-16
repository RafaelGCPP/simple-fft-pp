#include <numbers> // for std::numbers::pi
#include <array>
#include "constexpr_math.h"

namespace sfft
{
    template <typename T, int N>
    struct TwiddleGenerator
    {
        // Precompute twiddle factors at compile time using std::array
        static constexpr std::array<T, N / 2> twiddles = []()
        {
            std::array<T, N / 2> arr{};
            for (size_t i = 0; i < N / 2; ++i)
            {
                double angle = -2.0 * std::numbers::pi * i / N;
                arr[i] = T(sfft::math::cos(angle), sfft::math::sin(angle));
            }
            return arr;
        }();
    };

    template <typename BaseTwidGen, int Stride>
    struct StridedTwiddleGenerator
    {
        struct View
        {
            constexpr auto operator[](size_t i) const
            {
                return BaseTwidGen::twiddles[i * Stride];
            }
        };
        static constexpr View twiddles{};
    };

}