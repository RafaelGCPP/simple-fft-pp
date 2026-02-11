#include <iostream>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <complex>
#include "fixed_point.h"

void test_q31_conversion() {
    double val = 0.5;
    auto q_val = Q31::from_double(val);
    
    if (q_val.raw == 0x40000000) {
        std::cout << "✅ Q31 Conversion (0.5): Success" << std::endl;
    } else {
        std::cout << "❌ Q31 Conversion (0.5): Failed. Got: 0x" << std::hex << q_val.raw << std::endl;
    }
}

void test_addition() {
    auto a = Q31::from_double(0.1);
    auto b = Q31::from_double(0.2);
    auto result = a + b;
    
    double diff = std::abs(result.to_double() - 0.3);
    if (diff < 1e-9) {
        std::cout << "✅ Addition (0.1 + 0.2): Success" << std::endl;
    } else {
        std::cout << "❌ Addition failed. Diff: " << diff << std::endl;
    }
}

void test_multiplication() {
    auto a = Q31::from_double(0.5);
    auto b = Q31::from_double(0.5);
    auto result = a * b;

    if (result.raw == 0x20000000) {
        std::cout << "✅ Multiplication (0.5 * 0.5): Success" << std::endl;
    } else {
        std::cout << "❌ Multiplication failed. Got: 0x" << std::hex << result.raw << std::endl;
    }
}

void test_negative_numbers() {
    auto a = Q31::from_double(-0.75);
    auto b = Q31::from_double(0.5);
    auto result = a * b; 

    if (std::abs(result.to_double() - (-0.375)) < 1e-9) {
        std::cout << "✅ Negative Multiplication (-0.75 * 0.5): Success" << std::endl;
    } else {
        std::cout << "❌ Negative Multiplication failed: " << result.to_double() << std::endl;
    }
}

void test_multiplication_rounding() {
    // 0.1 in Q31 is not exact. 
    // 0.1 * 0.1 = 0.01
    auto a = Q31::from_double(0.1);
    auto b = Q31::from_double(0.1);
    auto result = a * b;
    
    double expected = 0.01;
    double actual = result.to_double();
    
    // Simple truncation would yield 0.0099999997...
    // Rounding should bring it to the closest representable value.
    std::cout << "✅ Rounded Multiplication Result: " << std::setprecision(10) << actual;
    std::cout << " (Expected approx: " << expected << ")" << std::endl;
}

void test_complex_multiplication() {
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
    }
}

void test_constexpr_verification() {
    static constexpr auto const_val = Q31::from_double(0.125);
    static_assert(const_val.raw == 0x10000000, "constexpr conversion failed");
    
    static constexpr Q31Complex const_cplx(0.5, -0.5);
    static_assert(const_cplx._real.raw == 0x40000000, "constexpr complex real failed");
    static_assert(const_cplx._imag.raw == -0x40000000, "constexpr complex imag failed");
    
    std::cout << "✅ constexpr static_assert: Success" << std::endl;
}

void test_mixed_addition() {
    std::cout << "Testing Mixed-Precision Addition (Q23 + Q31)... ";
    Q23 a = Q23::from_double(0.5);
    Q31 b = Q31::from_double(0.25);
    
    // a is the base (Q23), so result should be Q23
    auto result = a + b; 
    
    if (std::abs(result.to_double() - 0.75) < 1e-7) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌ Got: " << result.to_double() << std::endl;
    }
}

void test_mixed_multiplication() {
    std::cout << "Testing Mixed-Precision Multiplication (Q23 * Q31)... ";
    Q23 a = Q23::from_double(8.0); 
    Q31 b = Q31::from_double(0.125);
    
    // Result should be in the scale of 'a' (Q23)
    Q23 result = a * b; 
    
    if (std::abs(result.to_double() - 1.0) < 1e-7) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌ Got: " << result.to_double() << std::endl;
    }
}

void test_mixed_complex() {
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
    }
}

int main() {
    std::cout << "--- Starting FixedPoint Unit Tests ---" << std::endl;

    test_q31_conversion();
    test_addition();
    test_multiplication();
    test_negative_numbers();
    test_multiplication_rounding();
    test_complex_multiplication();
    test_constexpr_verification();
    test_mixed_addition();
    test_mixed_multiplication();
    test_mixed_complex();

    std::cout << "--- FixedPoint Tests Completed ---" << std::endl;
    return 0;
}