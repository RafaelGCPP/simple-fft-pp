#include <iostream>
#include <complex>
#include <cmath>
#include "fft_core.h"
#include "fixed_point.h"

using DoubleCplx = std::complex<double>;

// Helper for precision comparison
bool is_near(double a, double b, double epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

// 1. Floating-Point DIF Test
void test_float_dif() {
    std::cout << "Testing Floating-Point DIF Butterfly... ";
    DoubleCplx a(0.5, 0.2), b(0.3, -0.1), w(0.70710678, -0.70710678);
    
    DoubleCplx a_exp = a + b;
    DoubleCplx b_exp = (a - b) * w;

    DIF_Butterfly<DoubleCplx>::process(a, b, w);

    if (is_near(a.real(), a_exp.real()) && is_near(b.real(), b_exp.real())) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌" << std::endl;
    }
}

// 2. Floating-Point DIT Test
void test_float_dit() {
    std::cout << "Testing Floating-Point DIT Butterfly... ";
    DoubleCplx a(0.8, -0.2), b(0.1, 0.4), w(0.0, -1.0);
    
    DoubleCplx b_rot = b * w;
    DoubleCplx a_exp = a + b_rot;
    DoubleCplx b_exp = a - b_rot;

    DIT_Butterfly<DoubleCplx>::process(a, b, w);

    if (is_near(a.real(), a_exp.real()) && is_near(b.real(), b_exp.real())) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌" << std::endl;
    }
}

// 3. Fixed-Point DIF Test (Q31)
void test_fixed_dif() {
    std::cout << "Testing Fixed-Point DIF Butterfly (Q31)... ";
    Q31Complex a(0.5, 0.2), b(0.3, -0.1), w(0.70710678, -0.70710678);
    
    // Manual Reference
    DoubleCplx a_exp = DoubleCplx(0.5, 0.2) + DoubleCplx(0.3, -0.1);
    DoubleCplx b_exp = (DoubleCplx(0.5, 0.2) - DoubleCplx(0.3, -0.1)) * DoubleCplx(0.70710678, -0.70710678);

    DIF_Butterfly<Q31Complex>::process(a, b, w);

    if (is_near(a.to_double_real(), a_exp.real()) && is_near(b.to_double_real(), b_exp.real())) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌" << std::endl;
    }
}

// 4. Fixed-Point DIT Test (Q23+Q31)
void test_fixed_dit() {
    std::cout << "Testing Mixed-Precision DIT Butterfly (Data:Q23, Twiddle:Q31)... ";
    
    // Data has 8 bits of headroom (Q23)
    Q23Complex a(0.8, -0.2);
    Q23Complex b(0.1, 0.4);
    
    // Twiddle has max precision (Q31)
    Q31Complex w(0.0, -1.0); 
    
    // This will now use the template operator*
    DIT_Butterfly<Q23Complex>::process(a, b, w); 

    // Manual Reference for validation
    DoubleCplx b_rot = DoubleCplx(0.1, 0.4) * DoubleCplx(0.0, -1.0);
    DoubleCplx a_exp = DoubleCplx(0.8, -0.2) + b_rot;
    DoubleCplx b_exp = DoubleCplx(0.8, -0.2) - b_rot;

    if (is_near(a.to_double_real(), a_exp.real()) && is_near(b.to_double_real(), b_exp.real())) {
        std::cout << "✅" << std::endl;
    } else {
        std::cout << "❌" << std::endl;
    }
}

int main() {
    std::cout << "--- Starting Butterfly Logic Tests ---" << std::endl;
    test_float_dif();
    test_float_dit();
    test_fixed_dif();
    test_fixed_dit();
    std::cout << "--- All Butterfly Tests Completed ---" << std::endl;
    return 0;
}