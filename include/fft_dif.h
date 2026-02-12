#pragma once

#include <cstddef>
#include "fft_core.h"
#include "twiddles.h"

/**
 * @brief Compile-time unrolled Stage Runner
 */
template<typename T, typename TwidGen, size_t Stride, size_t NumGroups>
struct DIF_FFT_Stage {
    static constexpr void process(T* data) {
        // Use a strided version of the twiddle generator to automatically
        // handle the twiddle index calculation based on the current stage.
        // This simplifies butterfly processing by removing index arithmetic.
        using StridedGen = StridedTwiddleGenerator<TwidGen, NumGroups>;

        for (size_t g = 0; g < NumGroups; ++g) {
            size_t offset = g * Stride * 2;
            for (size_t i = 0; i < Stride; ++i) {
                auto w = StridedGen::get_twiddle(i);
                DIF_Butterfly<T>::process(data[offset + i], data[offset + i + Stride], w);
            }
        }

        // Tail Recursion: Move to the next stage
        // The 'if constexpr' ensures the recursion stops at compile time
        if constexpr (Stride > 1) {
            DIF_FFT_Stage<T, TwidGen, Stride / 2, NumGroups * 2>::process(data);
        }
    }
};


/**
 * @brief Forward FFT (DIF) with compile-time unrolling
 */
template<typename T, typename TwidGen>
class ForwardFFT_DIF {
public:
    static void process(T* data) {
        // Kick off the recursion starting at the largest stride
        DIF_FFT_Stage<T, TwidGen, TwidGen::N_Value / 2, 1>::process(data);
    }
};


