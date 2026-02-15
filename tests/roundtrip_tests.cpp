#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <complex>
#include "simple_fft.h"

using namespace sfft;

template<typename Complex, typename TwidComplex>
int test_fft_roundtrip()
{
    const size_t N = 64; // Larger N to test error accumulation
    const double TOLERANCE = 1e-4;

    std::cout << "--- FFT Roundtrip Test (N=" << N << ") ---" << std::endl;

    // 1. Type configuration

    auto fft = FFT<Complex, TwidComplex, N>();

    // 2. Create original signal (Sine + Cosine to test real and imag)
    std::vector<Complex> original(N);
    std::vector<Complex> buffer(N);

    for (size_t n = 0; n < N; ++n)
    {
        double val_r = 0.4 * std::sin(2.0 * M_PI * 2.0 * n / N); // 2 cycles
        double val_i = 0.3 * std::cos(2.0 * M_PI * 5.0 * n / N); // 5 cycles
        original[n] = Complex(val_r, val_i);
        buffer[n] = original[n];
    }

    // 3. Processing
    std::cout << "Step 1: Forward FFT (DIF)... " << std::flush;
    fft.process(buffer.data());
    std::cout << "Done." << std::endl;

    // Here the buffer is in Bit-Reversed order, but it doesn't matter!
    // The IFFT DIT expects exactly that.

    std::cout << "Step 2: Inverse FFT (DIT)... " << std::flush;
    fft.inverse(buffer.data());
    std::cout << "Done." << std::endl;

    // 4. Error Verification (MSE)
    double total_error = 0;
    bool failed = false;

    std::cout << "\nVerifying Results (First 8 samples):" << std::endl;
    std::cout << "n   | Original           | Recovered          | Diff" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    for (size_t n = 0; n < N; ++n)
    {
        double r_orig = (double)original[n].real();
        double i_orig = (double)original[n].imag();
        double r_rec = (double)buffer[n].real();
        double i_rec = (double)buffer[n].imag();
        
        double err_r = std::abs(r_orig - r_rec);
        double err_i = std::abs(i_orig - i_rec);
        total_error += (err_r * err_r) + (err_i * err_i);

        if (n < 8)
        {
            std::cout << std::setw(3) << n << " | "
                      << std::fixed << std::setprecision(4)
                      << r_orig << " + " << i_orig << "j | "
                      << r_rec << " + " << i_rec << "j | "
                      << std::scientific << std::setprecision(2) << std::max(err_r, err_i) << std::endl;
        }

        if (err_r > TOLERANCE || err_i > TOLERANCE)
        {
            failed = true;
        }
    }

    double mse = total_error / N;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Mean Squared Error: " << std::scientific << mse << std::endl;

    if (!failed)
    {
        std::cout << "\n✨ SUCCESS: Roundtrip verified within tolerance!" << std::endl;
    }
    else
    {
        std::cout << "\n❌ FAILURE: Signal distortion too high." << std::endl;
        return -1;
    }
    return 0;
}

int main()
{
    int retval = 0;
    retval |= test_fft_roundtrip<std::complex<double>, std::complex<double>>();
    retval |= test_fft_roundtrip<Q15Complex,Q15Complex>();
    retval |= test_fft_roundtrip<Q15Complex, Q31Complex>();
    return retval;
}