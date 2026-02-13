#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include "fixed_point.h"
#include "rfft.h"

using namespace sfft;

template <typename T, typename CplxT>
void test_known_complex_sequence()
{
    const size_t N = 16;
    std::cout << "\n--- Testing FFT Forward DIF (Known Sequence N=" << N << ") ---" << std::endl;
    T input[N] = {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
        -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0};

    T expected_real[N] = {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
        -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0};

    CplxT expected_fft[N] = {
        CplxT(0.00000, 0.00000),
        CplxT(9.00000, -45.24606),
        CplxT(-8.00000, 19.31371),
        CplxT(9.00000, -13.46945),
        CplxT(-8.00000, 8.00000),
        CplxT(9.00000, -6.01361),
        CplxT(-8.00000, 3.31371),
        CplxT(9.00000, -1.79021),
        CplxT(-8.00000, 0.00000),
        CplxT(9.00000, 1.79021),
        CplxT(-8.00000, -3.31371),
        CplxT(9.00000, 6.01361),
        CplxT(-8.00000, -8.00000),
        CplxT(9.00000, 13.46945),
        CplxT(-8.00000, -19.31371),
        CplxT(9.00000, 45.24606),
    };

    auto fft = RFFT<T, CplxT, TwiddleGenerator<CplxT, N>>();
    auto view = fft.process(input);

    std::cout << "Index | Real       | Imag       | Mag        | Expected Real | Expected Imag | Expected Mag | Abs Error" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    bool pass = true;
    for (size_t i = 0; i < N; ++i)
    {
        double r = (double)view[i].real();
        double img = (double)view[i].imag();
        double mag = std::sqrt(r * r + img * img);
        double exp_r = (double)expected_fft[i].real();
        double exp_img = (double)expected_fft[i].imag();
        double exp_mag = std::sqrt(exp_r * exp_r + exp_img * exp_img);
        double abs_error = std::sqrt((r - exp_r) * (r - exp_r) + (img - exp_img) * (img - exp_img));
        if (abs_error > 1e-4)
        {
            pass = false;
        }
        std::cout << std::setw(5) << i << " | "
                  << std::fixed << std::setprecision(5)
                  << std::setw(10) << r << " | "
                  << std::setw(10) << img << " | "
                  << std::setw(10) << mag << " | "
                  << std::setw(13) << exp_r << " | "
                  << std::setw(13) << exp_img << " | "
                  << std::setw(13) << exp_mag << " | "
                  << std::setw(9) << abs_error << std::endl;
    }

    auto ifft = IRFFT<T, CplxT, TwiddleGenerator<CplxT, N>>();
    T *output = ifft.process(input);
    std::cout << "\nIndex | Reconstructed Real | Expected Real | Abs Error" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    for (size_t i = 0; i < N; ++i)
    {
        double recon_r = (double)output[i];
        double exp_r = (double)expected_real[i];
        double abs_error = std::abs(recon_r - exp_r);
        if (abs_error > 1e-4)
        {
            pass = false;
        }
        std::cout << std::setw(5) << i << " | "
                  << std::fixed << std::setprecision(5)
                  << std::setw(18) << recon_r << " | "
                  << std::setw(13) << exp_r << " | "
                  << std::setw(9) << abs_error << std::endl;
    }

    if (pass)
    {
        std::cout << "✨ Test passed!" << std::endl;
    }
    else
    {
        std::cout << "❌ Test failed!" << std::endl;
    }
}
template <typename T, typename CplxT>
void test_rfft_filter_transform()
{
    const size_t N_REAL = 32;
    T buffer[N_REAL];

    // 1. Generate signal: Sine(2Hz) + Sine(12Hz)
    for (size_t i = 0; i < N_REAL; ++i)
    {
        double t = (double)i / N_REAL;
        double signal = 0.4 * std::sin(2.0 * M_PI * 2.0 * t) +
                        0.2 * std::sin(2.0 * M_PI * 12.0 * t);
        buffer[i] = T(signal);
    }

    // 2. Process RFFT
    using FullTwidGen = TwiddleGenerator<CplxT, N_REAL>;
    auto view = RFFT<T, CplxT, FullTwidGen>::process(buffer);

    // 3. Apply Transform (Low-Pass Filter)
    // If k > 5, zero out the bin.
    view.transform([](size_t k, CplxT val)
                   {
        if (k > 5) return CplxT(0.0, 0.0);
        return val; });

    // 4. Reconstruct via IRFFT
    IRFFT<T, CplxT, FullTwidGen>::process(buffer);

    bool pass = true;
    std::cout << std::endl;
    std::cout << "Reconstructed Signal (after filtering):" << std::endl;
    // 5. Verification
    std::cout << "Index | Filtered Output | Expected (2Hz only) | Error" << std::endl;
    for (size_t i = 0; i < N_REAL; ++i)
    {
        double t = (double)i / N_REAL;
        double expected = 0.4 * std::sin(2.0 * M_PI * 2.0 * t);
        double actual = (double)buffer[i];
        std::cout << std::setw(5) << i << " | "
                  << std::setw(15) << actual << " | "
                  << std::setw(18) << expected << " | "
                  << std::abs(actual - expected) << std::endl;
        if (std::abs(actual - expected) > 1e-4)
        {
            pass = false;
        }
    }

    if (pass)
    {
        std::cout << "✨ Test passed!" << std::endl;
    }
    else
    {
        std::cout << "❌ Test failed!" << std::endl;
    }
}
int main()
{
    test_known_complex_sequence<Q23, Q23Complex>();
    test_known_complex_sequence<double, std::complex<double>>();
    test_rfft_filter_transform<Q23, Q23Complex>();
    test_rfft_filter_transform<double, std::complex<double>>();
    return 0;
}