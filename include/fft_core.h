#pragma once

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