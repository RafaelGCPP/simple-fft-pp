#pragma once
#include <cstddef>
#include <algorithm> // for std::swap
#include "twiddles.h"

/**
 * @brief Radix-2 Decimation in Frequency (DIF) Butterfly
 * * Logic:
 * A_out = A + B
 * B_out = (A - B) * Twiddle
 * * Used in Forward FFT.
 */
template<typename T>
struct DIF_Butterfly {
    template<typename U>
    static constexpr void process(T& a, T& b, const U& twiddle) {
        T temp_a = a;
        T temp_b = b;

        a = temp_a + temp_b;
        T diff = temp_a - temp_b;
        
        // Complex multiplication rotation
        b = diff * twiddle;
    }
};

/**
 * @brief Radix-2 Decimation in Time (DIT) Butterfly
 * * Logic:
 * A_out = A + (B * Twiddle)
 * B_out = A - (B * Twiddle)
 * * Used in Inverse FFT.
 */
template<typename T>
struct DIT_Butterfly {
    template<typename U>
    static constexpr void process(T& a, T& b, const U& twiddle) {
        // Multiply first in DIT
        T rotated_b = b * twiddle;
        
        T temp_a = a;
        
        a = temp_a + rotated_b;
        b = temp_a - rotated_b;
    }
};



/**
 * @brief Compile-time bit reversal of an index
 */
template<size_t Bits>
constexpr size_t reverse_bits(size_t n) {
    size_t result = 0;
    for (size_t i = 0; i < Bits; ++i) {
        result = (result << 1) | (n & 1);
        n >>= 1;
    }
    return result;
}

/**
 * @brief Decorator that provides a linear view and an in-place commit method.
 */
template<typename T, size_t N>
class BitReversedView {
    T* const _data;
    static constexpr size_t _bits = std::countr_zero(N); 

public:
    explicit BitReversedView(T* data) : _data(data) {}

    // Accessor (The "Decorator" part)
    constexpr T& operator[](size_t index) {
        return _data[reverse_bits<_bits>(index)];
    }

    constexpr const T& operator[](size_t index) const {
        return _data[reverse_bits<_bits>(index)];
    }

    /**
     * @brief "Commits" the bit-reversal by physically reordering the data in RAM.
     * This makes the buffer linear so standard pointers/loops work.
     */
    void commit() {
        for (size_t i = 0; i < N; ++i) {
            size_t j = reverse_bits<_bits>(i);
            if (i < j) {
                std::swap(_data[i], _data[j]);
            }
        }
    }

    static constexpr size_t size() { return N; }
};

template<typename T>
constexpr T scale_in_half(const T& value) {
    return value * T(0.5);
}



/**
 * @brief Compile-time unrolled Stage Runner
 */
template<typename T, int N, typename TwidGen, typename Butterfly = DIF_Butterfly<T>, bool inverse = false>
struct FFT_Block {
    static constexpr void process(T* data) {    
        for (size_t i = 0; i < N/2 ; ++i) {
            auto w = TwidGen::get_twiddle(i);
            w = inverse ? conj(w) : w; // Conjugate for IFFT
            Butterfly::process(data[i], data[i + N/2], w);
            if (inverse) {
                data[i] = scale_in_half(data[i]);
                data[i + N/2] = scale_in_half(data[i + N/2]);
            }
        }
    }
};


template<typename T, typename TwidGen>
class FFT {
public:
    static auto process(T* data) {
        constexpr size_t N = TwidGen::N_Value;
        FFT_Block<T, N, TwidGen>::process(data);
        if constexpr (N > 2) {
            using StrideGen = StridedTwiddleGenerator<TwidGen, 2>;
            FFT<T, StrideGen>::process(data);
            FFT<T, StrideGen>::process(data + N/2);
        }
        return BitReversedView<T, N>(data); // Return a view for bit-reversal access
    }
};

template<typename T, typename TwidGen>
class IFFT {
public:
    static auto process(T* data) {
        constexpr size_t N = TwidGen::N_Value;
        if constexpr (N > 2) {
            using StrideGen = StridedTwiddleGenerator<TwidGen, 2>;
            IFFT<T, StrideGen>::process(data);
            IFFT<T, StrideGen>::process(data + N/2);
        }
        FFT_Block<T, N, TwidGen, DIT_Butterfly<T>, true>::process(data);
        return data; // Return raw pointer for linear access after IFFT
    }
};

