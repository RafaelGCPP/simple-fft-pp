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

    if (std::abs(result.to_double_real() - 0.5) < 1e-9 && 
        std::abs(result.to_double_imag() - 0.0) < 1e-9) {
        std::cout << "✅ Complex Multiplication (conjugates): Success" << std::endl;
    } else {
        std::cout << "❌ Complex Multiplication failed. Got: " 
                  << result.to_double_real() << " + " << result.to_double_imag() << "i" << std::endl;
    }
}

void test_constexpr_verification() {
    static constexpr auto const_val = Q31::from_double(0.125);
    static_assert(const_val.raw == 0x10000000, "constexpr conversion failed");
    
    static constexpr Q31Complex const_cplx(0.5, -0.5);
    static_assert(const_cplx.real.raw == 0x40000000, "constexpr complex real failed");
    
    std::cout << "✅ constexpr static_assert: Success" << std::endl;
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

    std::cout << "--- FixedPoint Tests Completed ---" << std::endl;
    return 0;
}