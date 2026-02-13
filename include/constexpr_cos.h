#pragma once

/**
 * @brief Constexpr-safe cosine using a Maclaurin series with range reduction.
 *
 * Since std::cos is not constexpr until C++26, this provides a portable
 * compile-time alternative. The input is first reduced to [0, pi/2] via
 * quadrant symmetry, then a truncated Maclaurin series is evaluated.
 * 12 terms give approximately 18 digits of precision for double.
 */
namespace sfft {

namespace detail {

constexpr double constexpr_abs(double x) {
    return x < 0.0 ? -x : x;
}

constexpr double constexpr_fmod(double x, double y) {
    return x - static_cast<long long>(x / y) * y;
}

/// Maclaurin series for cos(x), accurate for small |x|.
/// cos(x) = sum_{n=0}^{N} (-1)^n * x^{2n} / (2n)!
constexpr double cos_series(double x) {
    const double x2 = x * x;
    double term = 1.0;
    double sum  = 1.0;
    for (int n = 1; n <= 12; ++n) {
        term *= -x2 / static_cast<double>((2 * n - 1) * (2 * n));
        sum += term;
    }
    return sum;
}

} // namespace detail

constexpr double constexpr_pi = 3.14159265358979323846;

/**
 * @brief Compute cos(x) at compile time.
 * @param x Angle in radians (any value).
 * @return cos(x) with ~18 digits of double precision.
 */
constexpr double constexpr_cos(double x) {
    // 1. Make x positive (cos is even)
    x = detail::constexpr_abs(x);

    // 2. Reduce to [0, 2*pi)
    x = detail::constexpr_fmod(x, 2.0 * constexpr_pi);

    // 3. Quadrant reduction to [0, pi/2]
    //    Q0: [0, pi/2)       -> cos_series(x)
    //    Q1: [pi/2, pi)      -> -cos_series(pi - x)
    //    Q2: [pi, 3*pi/2)    -> -cos_series(x - pi)
    //    Q3: [3*pi/2, 2*pi)  ->  cos_series(2*pi - x)
    constexpr double half_pi      = constexpr_pi / 2.0;
    constexpr double three_half_pi = 3.0 * constexpr_pi / 2.0;

    if (x <= half_pi) {
        return detail::cos_series(x);
    } else if (x <= constexpr_pi) {
        return -detail::cos_series(constexpr_pi - x);
    } else if (x <= three_half_pi) {
        return -detail::cos_series(x - constexpr_pi);
    } else {
        return detail::cos_series(2.0 * constexpr_pi - x);
    }
}

} // namespace sfft