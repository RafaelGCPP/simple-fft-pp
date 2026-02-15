#pragma once

#include <limits>

namespace sfft {
namespace math {

constexpr double pi = 3.14159265358979323846;

// Constexpr implementation of power function
constexpr double pow(double base, int exp) {
    double res = 1.0;
    for (int i = 0; i < exp; ++i) res *= base;
    return res;
}

// Constexpr implementation of factorial
constexpr double factorial(int n) {
    double res = 1.0;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

// Normalize angle to range [-pi, pi]
constexpr double normalize_angle(double x) {
    // Basic range reduction
    if (x > pi || x < -pi) {
        long long k = static_cast<long long>(x / (2 * pi));
        x -= k * 2 * pi;
    }
    while (x > pi) x -= 2 * pi;
    while (x < -pi) x += 2 * pi;
    return x;
}

// Taylor sine series
constexpr double sin(double x) {
    x = normalize_angle(x);
    double res = x;
    double term = x;
    double x2 = x * x;
    
    // 15 terms should be enough for double precision
    for (int i = 1; i < 15; ++i) {
        term *= -x2 / ((2 * i) * (2 * i + 1));
        res += term;
    }
    return res;
}

// Taylor series for cos(x)
// x should be in range [-pi, pi] for better convergence
constexpr double cos(double x) {
    x = normalize_angle(x);
    double res = 1.0;
    double term = 1.0;
    double x2 = x * x;
    
    // 15 terms should be enough for double precision
    for (int i = 1; i < 15; ++i) {
        term *= -x2 / ((2 * i - 1) * (2 * i));
        res += term;
    }
    return res;
}

} // namespace math
} // namespace sfft
