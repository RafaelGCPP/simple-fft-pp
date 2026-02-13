#include <iostream>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <complex>
#include "fixed_point.h"

using namespace sfft;

int test_q31_conversion() {
    double val = 0.5;
    auto q_val = Q31(val);
    
    if (q_val.raw == 0x40000000) {
        std::cout << "✅ Q31 Conversion (0.5): Success" << std::endl;
    } else {
        std::cout << "❌ Q31 Conversion (0.5): Failed. Got: 0x" << std::hex << q_val.raw << std::endl;
        return -1;
    }
    return 0;
}

int test_addition() {
    auto a = Q31(0.1);
    auto b = Q31(0.2);
    auto result = a + b;
    
    double diff = std::abs((double)result - 0.3);
    if (diff < 1e-9) {
        std::cout << "✅ Addition (0.1 + 0.2): Success" << std::endl;
    } else {
        std::cout << "❌ Addition failed. Diff: " << diff << std::endl;
        return -1;
    }
    return 0;
}

int test_multiplication() {
    auto a = Q31(0.5);
    auto b = Q31(0.5);
    auto result = a * b;

    if (result.raw == 0x20000000) {
        std::cout << "✅ Multiplication (0.5 * 0.5): Success" << std::endl;
    } else {
        std::cout << "❌ Multiplication failed. Got: 0x" << std::hex << result.raw << std::endl;
        return -1;
    }
    return 0;
}

int test_negative_numbers() {
    auto a = Q31(-0.75);
    auto b = Q31(0.5);
    auto result = a * b; 

    if (std::abs((double)result - (-0.375)) < 1e-9) {
        std::cout << "✅ Negative Multiplication (-0.75 * 0.5): Success" << std::endl;
    } else {
        std::cout << "❌ Negative Multiplication failed: " << (double)result << std::endl;
        return -1;
    }
    return 0;
}

int test_multiplication_rounding() {
    // 0.1 in Q31 is not exact. 
    // 0.1 * 0.1 = 0.01
    auto a = Q31(0.1);
    auto b = Q31(0.1);
    auto result = a * b;
    
    double expected = 0.01;
    double actual = (double)result;
    
    // Simple truncation would yield 0.0099999997...
    // Rounding should bring it to the closest representable value.
    std::cout << "✅ Rounded Multiplication Result: " << std::setprecision(10) << actual;
    std::cout << " (Expected approx: " << expected << ")" << std::endl;
    return 0;
}

int test_complex_multiplication() {
    // (0.5 + 0.5i) * (0.5 - 0.5i) = 0.25 - (-0.25) + ( -0.25 + 0.25)i = 0.5 + 0i
    Q31Complex a(0.5, 0.5);
    Q31Complex b(0.5, -0.5);
    auto result = a * b;

    if (std::abs((double)result.real() - 0.5) < 1e-9 && 
        std::abs((double)result.imag() - 0.0) < 1e-9) {
        std::cout << "✅ Complex Multiplication (conjugates): Success" << std::endl;
    } else {
        std::cout << "❌ Complex Multiplication failed. Got: " 
                  << (double)result.real() << " + " << (double)result.imag() << "i" << std::endl;
        return -1;
    }
    return 0;
}

int test_constexpr_verification() {
    static constexpr auto const_val = Q31(0.125);
    static_assert(const_val.raw == 0x10000000, "constexpr conversion failed");
    
    static constexpr Q31Complex const_cplx(0.5, -0.5);
    static_assert(const_cplx._real.raw == 0x40000000, "constexpr complex real failed");
    static_assert(const_cplx._imag.raw == -0x40000000, "constexpr complex imag failed");
    
    std::cout << "✅ constexpr static_assert: Success" << std::endl;
    return 0;
}

int test_mixed_addition() {
    std::cout << "Testing Mixed-Precision Addition (Q23 + Q31)... ";
    Q23 a = Q23(0.5);
    Q31 b = Q31(0.25);
    
    // a is the base (Q23), so result should be Q23
    auto result = a + b; 
    
    if (std::abs((double)result - 0.75) < 1e-7) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌ Got: " << (double)result << std::endl;
        return -1;
    }
    return 0;
}

int test_mixed_multiplication() {
    std::cout << "Testing Mixed-Precision Multiplication (Q23 * Q31)... ";
    Q23 a = Q23(8.0); 
    Q31 b = Q31(0.125);
    
    // Result should be in the scale of 'a' (Q23)
    Q23 result = a * b; 
    
    if (std::abs((double)result - 1.0) < 1e-7) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌ Got: " << (double)result << std::endl;
        return -1;
    }
    return 0;
}

int test_mixed_complex() {
    std::cout << "Testing Mixed-Precision Complex Multiply (Q23Cplx * Q31Cplx)... ";
    // Common case for our FFT: Signal in Q23, Twiddle in Q31
    Q23Complex signal(0.5, 0.0);
    Q31Complex twiddle(0.0, -1.0); // -j
    
    auto result = signal * twiddle; // Should result in (0, -0.5) in Q23 scale

    if (std::abs((double)result.real() - 0.0) < 1e-7 && 
        std::abs((double)result.imag() - (-0.5)) < 1e-7) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌ Got: " << (double)result.real() << " + " << (double)result.imag() << "i" << std::endl;
        return -1;
    }
    return 0;
}

int main() {
    std::cout << "--- Starting FixedPoint Unit Tests ---" << std::endl;
    int retval=0;

    retval |= test_q31_conversion();
    retval |= test_addition();
    retval |= test_multiplication();
    retval |= test_negative_numbers();
    retval |= test_multiplication_rounding();
    retval |= test_complex_multiplication();
    retval |= test_constexpr_verification();
    retval |= test_mixed_addition();
    retval |= test_mixed_multiplication();
    retval |= test_mixed_complex();

    std::cout << "--- FixedPoint Tests Completed ---" << std::endl;
    return retval;
}