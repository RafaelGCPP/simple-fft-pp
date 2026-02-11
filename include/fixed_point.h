#pragma once

#include <cstdint>
#include <cmath>

/**
 * @brief A generic Fixed Point scalar type.
 * @tparam FractionalBits Number of bits used for the fractional part (e.g., 31 for Q31).
 */
template<int FractionalBits>
struct FixedPoint {
    int32_t raw;

    // Default constructor
    constexpr FixedPoint() : raw(0) {}
    
    // Explicit constructor from raw bits
    constexpr explicit FixedPoint(int32_t r) : raw(r) {}

    // Factory from double (calculated at compile-time if constexpr)
    static constexpr FixedPoint from_double(double d) {
    // Calculate the maximum possible positive value for this bit depth
    // For Q31 (31 fractional bits), this is 0x7FFFFFFF
    constexpr int32_t kMaxPos = std::numeric_limits<int32_t>::max();
    constexpr int32_t kMaxNeg = std::numeric_limits<int32_t>::min();

    double scaled = d * (1LL << FractionalBits);
    
    if (scaled >= static_cast<double>(kMaxPos)) return FixedPoint(kMaxPos);
    if (scaled <= static_cast<double>(kMaxNeg)) return FixedPoint(kMaxNeg);
    
    return FixedPoint(static_cast<int32_t>(scaled));
    }

    // Converters
    constexpr double to_double() const {
        return static_cast<double>(raw) / (1LL << FractionalBits);
    }

    // --- Arithmetic Operators ---

    constexpr FixedPoint operator+(const FixedPoint& other) const {
        return FixedPoint(raw + other.raw);
    }

    constexpr FixedPoint operator-(const FixedPoint& other) const {
        return FixedPoint(raw - other.raw);
    }

    constexpr FixedPoint operator-() const {
        return FixedPoint(-raw);
    }

    /**
    * @brief Multiplication with Rounding to Nearest.
    */
    constexpr FixedPoint operator*(const FixedPoint& other) const {
        int64_t product = static_cast<int64_t>(raw) * other.raw;
        
        // Add rounding bias: 1 << (FractionalBits - 1)
        // This represents 0.5 in the discarded fractional part.
        int64_t rounded_product = product + (1LL << (FractionalBits - 1));
        
        return FixedPoint(static_cast<int32_t>(rounded_product >> FractionalBits));
    }
    // Support for mixed-precision multiplication could be added later
};

// Common aliases
using Q31 = FixedPoint<31>;
using Q23 = FixedPoint<23>;

/**
 * @brief A generic complex number template using fixed-point types.
 * @tparam T A fixed-point type (e.g., Q31, Q23).
 */
template<typename T>
struct FixedComplex {
    T real;
    T imag;

    constexpr FixedComplex() : real(), imag() {}
    constexpr FixedComplex(T r, T i) : real(r), imag(i) {}
    
    // Explicitly handle construction to avoid "explicit" constructor conflicts
    constexpr FixedComplex(double r, double i) 
        : real(T::from_double(r)), imag(T::from_double(i)) {}

    // --- Complex Arithmetic ---

    constexpr FixedComplex operator+(const FixedComplex& other) const {
        return {real + other.real, imag + other.imag};
    }

    constexpr FixedComplex operator-(const FixedComplex& other) const {
        return {real - other.real, imag - other.imag};
    }

    /**
     * @brief Complex Multiplication (The Twiddle Rotation)
     * (a+bi)*(c+di) = (ac-bd) + (ad+bc)i
     */
    constexpr FixedComplex operator*(const FixedComplex& other) const {
        return {
            (real * other.real) - (imag * other.imag),
            (real * other.imag) + (imag * other.real)
        };
    }

    // Helper for debugging/testing
    constexpr double to_double_real() const { return real.to_double(); }
    constexpr double to_double_imag() const { return imag.to_double(); }
};

// Common aliases
using Q31Complex = FixedComplex<Q31>;
using Q23Complex = FixedComplex<Q23>;
