#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <vector>
#include "simple_fft.h"

using namespace sfft;

// Local helper for testing purposes only
double get_abs_diff(Q31Complex actual, std::complex<double> expected) {
    double real_diff = (double) actual.real() - expected.real();
    double imag_diff = (double) actual.imag() - expected.imag();
    return std::sqrt(real_diff * real_diff + imag_diff * imag_diff);
}

// Helper to print failures clearly
void print_failure(int k, std::complex<double> expected, std::complex<double> actual) {
    std::cout << "❌ Error at index k=" << k << "\n"
              << "   Expected: " << expected.real() << " + " << expected.imag() << "i\n"
              << "   Got:      " << (double) actual.real() << " + " << (double) actual.imag() << "i\n"
              << "   Diff:     " << std::abs(expected - actual) << "\n\n";
}

// Helper to print Q31 failures
void print_q31_failure(int k, std::complex<double> expected, Q31Complex actual) {
    std::cout << "❌ Error at index k=" << k << "\n"
              << "   Expected: " << expected.real() << " + " << expected.imag() << "i\n"
              << "   Got:      " << (double) actual.real() << " + " <<  (double) actual.imag() << "i\n"
              << "   Diff:     " << get_abs_diff(actual, expected) << "\n\n";
}

void test_twiddle_generation() {
    const int N = 256;
    const double EPSILON = 1e-12; 
    bool all_passed = true;

    using TestGen = TwiddleGenerator<std::complex<double>, N>;

    std::cout << "--- Starting Twiddle Test (N=" << N << ") ---\n";

    for (int k = 0; k <= N / 2; ++k) {
        double angle = -2.0 * M_PI * k / N;
        std::complex<double> expected = std::polar(1.0, angle);
        std::complex<double> actual = TestGen::get_twiddle(k);

        if (std::abs(expected - actual) > EPSILON) {
            print_failure(k, expected, actual);
            all_passed = false;
        } else {
            std::cout << "✅ k=" << std::setw(3) << k << " OK ";
            if ((k + 1) % 8 == 0) std::cout << std::endl;
        }
    }

    if (all_passed) {
        std::cout << "\n✨ All twiddles are mathematically correct!\n";
    } else {
        std::cout << "\n⚠️ The test failed at some points.\n";
    }
}

void test_twiddle_generation_q31() {
    const int N = 256;
    const double EPSILON = 1e-6; 
    bool all_passed = true;

    using TestGen = TwiddleGenerator<Q31Complex, N>;

    std::cout << "\n--- Starting Q31 Fixed-Point Twiddle Test (N=" << N << ") ---\n";

    for (int k = 0; k <= N / 2; ++k) {
        double angle = -2.0 * M_PI * k / N;
        std::complex<double> expected = std::polar(1.0, angle);
        Q31Complex actual = TestGen::get_twiddle(k);

        double diff = get_abs_diff(actual, expected);

        if (diff > EPSILON) {
            print_q31_failure(k, expected, actual);
            all_passed = false;
        } else {
            std::cout << "✅ k=" << std::setw(3) << k << " OK ";
            if ((k + 1) % 8 == 0) std::cout << std::endl;
        }
    }

    if (all_passed) {
        std::cout << "\n✨ All Q31 twiddles are within tolerance!\n";
    } else {
        std::cout << "\n⚠️ The Q31 test failed at some points.\n";
    }
}

int main() {
    test_twiddle_generation();
    test_twiddle_generation_q31();
    return 0;
}