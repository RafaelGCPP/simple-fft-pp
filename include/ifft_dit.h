#pragma once

#include "fft_core.h"
#include "twiddles.h"
#include <cstddef>


/**
 * @brief Generic function to calculate half of a value.
 * Default version for types that support multiplication by double (like std::complex).
 */
template<typename T>
constexpr T scale_in_half(const T& value) {
    return value * T(0.5);
}


template <typename T, typename TwidGen, size_t Stride, size_t NumGroups>
struct IFFT_Stage
{
    static constexpr void process(T *data)
    {
        // Use a strided version of the twiddle generator to automatically
        // handle the twiddle index calculation based on the current stage.
        using StridedGen = StridedTwiddleGenerator<TwidGen, NumGroups>;

        for (size_t g = 0; g < NumGroups; ++g)
        {
            size_t offset = g * Stride * 2;
            for (size_t i = 0; i < Stride; ++i)
            {
                auto w = conj(StridedGen::get_twiddle(i));

                DIT_Butterfly<T>::process(
                    data[offset + i],
                    data[offset + i + Stride],
                    w);

                data[offset + i] = scale_in_half(data[offset + i]);
                data[offset + i + Stride] = scale_in_half(data[offset + i + Stride]);
            }
        }

        if constexpr (Stride * 2 < TwidGen::N_Value)
        {
            IFFT_Stage<T, TwidGen, Stride * 2, NumGroups / 2>::process(data);
        }
    }
};

template <typename T, typename TwidGen>
class InverseFFT_DIT
{
public:
    static void process(T *data)
    {
        IFFT_Stage<T, TwidGen, 1, TwidGen::N_Value / 2>::process(data);
    }
};

