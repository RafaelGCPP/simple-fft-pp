#pragma once

#include <cstdint>
#include <cmath>

/**
 * @brief A generic Fixed Point scalar type.
 * @tparam B Number of bits used for the fractional part (e.g., 31 for Q31).
 */
template <int B>
struct FixedPoint
{
    static constexpr int FractionalBits = B; // Adicione isso!
    int32_t raw;

    // Default constructor
    constexpr FixedPoint() : raw(0) {}

    // Explicit constructor from raw bits
    constexpr explicit FixedPoint(int32_t r) : raw(r) {}

    // Factory from double (calculated at compile-time if constexpr)
    static constexpr FixedPoint from_double(double d)
    {
        // Calculate the maximum possible positive value for this bit depth
        // For Q31 (31 fractional bits), this is 0x7FFFFFFF
        constexpr int32_t kMaxPos = std::numeric_limits<int32_t>::max();
        constexpr int32_t kMaxNeg = std::numeric_limits<int32_t>::min();

        double scaled = d * (1LL << FractionalBits);

        if (scaled >= static_cast<double>(kMaxPos))
            return FixedPoint(kMaxPos);
        if (scaled <= static_cast<double>(kMaxNeg))
            return FixedPoint(kMaxNeg);

        return FixedPoint(static_cast<int32_t>(scaled));
    }

    // Converters
    constexpr double to_double() const
    {
        return static_cast<double>(raw) / (1LL << FractionalBits);
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

    // Operador de conversão explícito para double
    explicit constexpr operator double() const {
        return static_cast<double>(raw) / (1LL << FractionalBits);
    }

    // Operador de conversão explícito para float
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
        : _real(T::from_double(r)), _imag(T::from_double(i)) {}

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
        // 1. Fazemos as 4 multiplicações em 64 bits puras (sem shift ainda)
        int64_t r1 = static_cast<int64_t>(this->_real.raw) * other.real().raw;
        int64_t r2 = static_cast<int64_t>(this->_imag.raw) * other.imag().raw;
        int64_t i1 = static_cast<int64_t>(this->_real.raw) * other.imag().raw;
        int64_t i2 = static_cast<int64_t>(this->_imag.raw) * other.real().raw;

        // 2. Soma/Subtrai no acumulador de 64 bits (precisão total)
        int64_t real_acc = r1 - r2;
        int64_t imag_acc = i1 + i2;

        // 3. Arredonda e faz o shift único pelos bits do Twiddle (OtherB)
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

    // Sobrecargas de função membro para imitar std::complex
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
 * @brief Sobrecarga específica para o nosso FixedComplex.
 * Aqui usamos o shift direto para máxima performance no RP2350.
 */
template<int B>
constexpr FixedComplex<FixedPoint<B>> scale_in_half(const FixedComplex<FixedPoint<B>>& value) {
    // Usamos o operador de shift que definimos na FixedPoint
    std::cout << "Scaling down by half: " << (double)value.real() << " + " << (double)value.imag() << "i -> ";
    auto result = FixedComplex<FixedPoint<B>>(
        FixedPoint<B>(value._real.raw + 1 >> 1),
        FixedPoint<B>(value._imag.raw + 1 >> 1)
    );
    std::cout << (double)result.real() << " + " << (double)result.imag() << "i" << std::endl;
    return result;
}