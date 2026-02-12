#pragma once

#include <cstdint>
#include <array>
#include <cmath>

/**
 * @brief Generic Twiddle generator using quadrant symmetry.
 * @tparam T The complex type (must accept construction by Real and Imaginary parts).
 * @tparam N FFT size.
 */
template<typename T, size_t N>
class TwiddleGenerator {
public:
    static constexpr size_t N_Value = N;
    static constexpr size_t Q = N / 4;

    // Base table (always stored in double or Q31 for maximum internal precision)
    struct InternalTable {
        std::array<double, Q + 1> cos_values;
        constexpr InternalTable() : cos_values{} {
            for (size_t i = 0; i <= Q; ++i) {
                // We use C++20 compile-time function if available,
                // or just constant calculations in modern compilers.
                cos_values[i] = std::cos((2.0 * M_PI * i) / N);
            }
        }
    };

    static constexpr InternalTable base_table = InternalTable();

    /**
     * @brief Factory method to create type T from double components.
     * This is the "glue" that allows using any type.
     */
    static constexpr T create_complex(double r, double i) {
        // If T is std::complex, this works.
        // If it's your FixedCplx, it will need a constructor that accepts double
        // or you can specialize this function.
        return T(r, i);
    }

    static constexpr T get_twiddle(size_t t) {
        bool is_2nd = (t > Q);
        
        size_t offset = is_2nd ? (2 * Q - 2 * t) : 0;
        size_t idx_r = t + offset;
        size_t idx_i = Q - t - offset;

        double r = base_table.cos_values[idx_r];
        double i = base_table.cos_values[idx_i];

        // Symmetry application:
        // 1st Q: (cos(t), -sin(t)) -> (table[idx_r], -table[idx_i])
        // 2nd Q: (-sin(t-Q), -cos(t-Q)) -> (-table[idx_r], -table[idx_i])
        
        double final_r = is_2nd ? -r : r;
        double final_i = -i; // The sine is always positive in 1st and 2nd Q, so -sin is always negative

        return create_complex(final_r, final_i);
    }

};

/**
 * @brief Visitor that applies a stride to twiddle factor access.
 * 
 * This wrapper allows a TwiddleGenerator to be used with different strides,
 * which is essential for enabling recursive FFT decomposition where each level
 * skips twiddle factors according to stage requirements.
 * 
 * Example usage:
 *   using BaseGen = TwiddleGenerator<Q31Complex, 16>;
 *   using StridedGen = StridedTwiddleGenerator<BaseGen, 2>;
 *   // Now StridedGen::get_twiddle(k) returns BaseGen::get_twiddle(2*k)
 * 
 * @tparam WrappedGen The base TwiddleGenerator type
 * @tparam Stride The stride factor to apply
 */
template<typename WrappedGen, size_t Stride>
class StridedTwiddleGenerator {
public:
    static constexpr size_t N_Value = WrappedGen::N_Value;
    static constexpr size_t S = Stride;

    /**
     * @brief Get twiddle factor at index t*Stride
     */
    static constexpr auto get_twiddle(size_t t) {
        return WrappedGen::get_twiddle(t * Stride);
    }
};

