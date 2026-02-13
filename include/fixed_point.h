#pragma once

#include <cstdint>
#include <cmath>

namespace sfft {

/**
 * @brief A generic Fixed Point scalar type.
 * @tparam B Number of bits used for the fractional part (e.g., 31 for Q31).
 */
template <int B>
struct FixedPoint
{
    static constexpr int FractionalBits = B;
    int32_t raw;

    // Default constructor
    constexpr FixedPoint() : raw(0) {}

    // Explicit constructor from raw bits
    constexpr explicit FixedPoint(int32_t r) : raw(r) {}

    // Constructor from double (replaces from_double)
    constexpr FixedPoint(double d) : raw(0)
    {
        constexpr int32_t kMaxPos = std::numeric_limits<int32_t>::max();
        constexpr int32_t kMaxNeg = std::numeric_limits<int32_t>::min();

        double scaled = d * (1LL << FractionalBits);

        if (scaled >= static_cast<double>(kMaxPos))
            raw = kMaxPos;
        else if (scaled <= static_cast<double>(kMaxNeg))
            raw = kMaxNeg;
        else
            raw = static_cast<int32_t>(scaled);
    }

    // --- Arithmetic Operators ---

    /**
     * @brief Mixed-precision Addition.
     * Aligns the 'other' scale to 'this' scale.
     */
    template <int OtherFracBits>
    constexpr FixedPoint operator+(const FixedPoint<OtherFracBits> &other) const
    {
        int32_t aligned_other;
        constexpr int shift = OtherFracBits - FractionalBits;

        if constexpr (shift > 0)
        {
            aligned_other = other.raw >> shift;
        }
        else if constexpr (shift < 0)
        {
            aligned_other = other.raw << (-shift);
        }
        else
        {
            aligned_other = other.raw;
        }

        return FixedPoint(raw + aligned_other);
    }

    /**
     * @brief Mixed-precision Subtraction.
     */
    template <int OtherFracBits>
    constexpr FixedPoint operator-(const FixedPoint<OtherFracBits> &other) const
    {
        int32_t aligned_other;
        constexpr int shift = OtherFracBits - FractionalBits;

        if constexpr (shift > 0)
        {
            aligned_other = other.raw >> shift;
        }
        else if constexpr (shift < 0)
        {
            aligned_other = other.raw << (-shift);
        }
        else
        {
            aligned_other = other.raw;
        }

        return FixedPoint(raw - aligned_other);
    }

    constexpr FixedPoint operator-() const
    {
        if (raw == std::numeric_limits<int32_t>::min()) {
            // If we negate -1, which is -2147483648 in two's complement,
            // we would overflow. So we clamp to max positive value.
            return FixedPoint(std::numeric_limits<int32_t>::max());
        }
        return FixedPoint(-raw);
    }

    /**
     * @brief Mixed-precision Multiplication.
     * Multiplies this (FractionalBits) by other (OtherFracBits).
     * Returns a result in the current FractionalBits.
     */
    template <int OtherFracBits>
    constexpr FixedPoint operator*(const FixedPoint<OtherFracBits> &other) const
    {
        int64_t product = static_cast<int64_t>(raw) * other.raw;
        // We must shift by OtherFracBits to return to our scale
        int64_t rounded = product + (1LL << (OtherFracBits - 1));
        return FixedPoint(static_cast<int32_t>(rounded >> OtherFracBits));
    }
    // Support for mixed-precision multiplication could be added later

    // Explicit conversion operator to double
    explicit constexpr operator double() const {
        return static_cast<double>(raw) / (1LL << FractionalBits);
    }

    // Explicit conversion operator to float
    explicit constexpr operator float() const {
        return static_cast<float>(raw) / (1LL << FractionalBits);
    }
};

// Common aliases
using Q31 = FixedPoint<31>;
using Q23 = FixedPoint<23>;
using Q15 = FixedPoint<15>;

/**
 * @brief A generic complex number template using fixed-point types.
 * @tparam T A fixed-point type (e.g., Q31, Q23).
 */
template <typename T>
struct FixedComplex
{
    T _real;
    T _imag;

    constexpr FixedComplex() : _real(), _imag() {}
    constexpr FixedComplex(T r, T i) : _real(r), _imag(i) {}

    // Explicitly handle construction to avoid "explicit" constructor conflicts
    constexpr FixedComplex(double r, double i)
        : _real(r), _imag(i) {}

    // --- Complex Arithmetic ---
    template <typename U>
    constexpr auto operator+(const FixedComplex<U> &other) const
    {
        // This relies on the FixedPoint template operator+ we just wrote
        return FixedComplex<T>{_real + other._real, _imag + other._imag};
    }

    template <typename U>
    constexpr auto operator-(const FixedComplex<U> &other) const
    {
        return FixedComplex<T>{_real - other._real, _imag - other._imag};
    }

    template<int OtherB>
    constexpr auto operator*(const FixedComplex<FixedPoint<OtherB>>& other) const {
        // 1. Perform 4 multiplications in pure 64 bits (no shift yet)
        int64_t r1 = static_cast<int64_t>(this->_real.raw) * other.real().raw;
        int64_t r2 = static_cast<int64_t>(this->_imag.raw) * other.imag().raw;
        int64_t i1 = static_cast<int64_t>(this->_real.raw) * other.imag().raw;
        int64_t i2 = static_cast<int64_t>(this->_imag.raw) * other.real().raw;

        // 2. Sum/Subtract in 64-bit accumulator (full precision)
        int64_t real_acc = r1 - r2;
        int64_t imag_acc = i1 + i2;

        // 3. Round and apply single shift by Twiddle bits (OtherB)
        real_acc += (1LL << (OtherB - 1));
        imag_acc += (1LL << (OtherB - 1));

        return FixedComplex<T>(
            FixedPoint<T::FractionalBits>(static_cast<int32_t>(real_acc >> OtherB)),
            FixedPoint<T::FractionalBits>(static_cast<int32_t>(imag_acc >> OtherB))
        );
    }


    constexpr FixedComplex conj() const
    {
        return {_real, -_imag};
    }

    friend constexpr FixedComplex conj(const FixedComplex &c)
    {
        return c.conj();
    }

    // Member function overloads to mimic std::complex
    constexpr T& real() { return _real; }
    constexpr const T& real() const { return _real; }
    constexpr T& imag() { return _imag; }
    constexpr const T& imag() const { return _imag; }
};

template<typename T>
constexpr auto real(const FixedComplex<T>& c) { return c.real(); }

template<typename T>
constexpr auto imag(const FixedComplex<T>& c) { return c.imag(); }

// Common aliases
using Q31Complex = FixedComplex<Q31>;
using Q23Complex = FixedComplex<Q23>;
using Q15Complex = FixedComplex<Q15>;




/**
 * @brief Specific overload for our FixedComplex.
 * Here we use direct bit-shift for maximum performance on RP2350.
 */
template<int B>
constexpr FixedComplex<FixedPoint<B>> scale_in_half(const FixedComplex<FixedPoint<B>>& value) {
    // Use the shift operator we defined in FixedPoint
    return FixedComplex<FixedPoint<B>>(
        FixedPoint<B>((value._real.raw + 1) >> 1),
        FixedPoint<B>((value._imag.raw + 1) >> 1)
    );
}

} // namespace sfft