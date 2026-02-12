#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include "fixed_point.h"
#include "rfft.h"

template <typename T, typename CplxT>
void test_known_complex_sequence()
{
    const size_t N = 16;
    std::cout << "\n--- Testing FFT Forward DIF (Known Sequence N=" << N << ") ---" << std::endl;
    T input[N] = {
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
    for (size_t i = 0; i < N ; ++i)
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
    if (pass)
    {
        std::cout << "✨ Test passed!" << std::endl;
    }
    else
    {
        std::cout << "❌ Test failed!" << std::endl;
    }
}

// template <typename T>
// void test_ifft_known_complex_sequence()
// {
//     const size_t N = 8;
//     std::cout << "\n--- Testing IFFT (Known Sequence N=" << N << ") ---" << std::endl;
//     T expected_result[N] = {
//         T(1.0, 2.0),
//         T(3.0, 4.0),
//         T(5.0, 6.0),
//         T(7.0, 8.0),
//         T(-8.0, -7.0),
//         T(-6.0, -5.0),
//         T(-4.0, -3.0),
//         T(-2.0, -1.0)
//     };

//     T input[N] = {
//         T(-4.0, 4.0),
//         T(30.72792, -12.72792),
//         T(-16.0, 0.0),
//         T(12.72792, 5.27208),
//         T(-8.0, -8.0),
//         T(5.27208, 12.72792),
//         T(0.0, -16.0),
//         T(-12.72792, 30.72792)
//     };

//     auto twid_gen = TwiddleGenerator<T, N>();
//     auto fft = IFFT<T, decltype(twid_gen)>();
//     auto view = BitReversedView<T, N>(input);
//     view.commit(); // Ensure data is bitreversed for IFFT

//     fft.process(input);
//     std::cout << "Index | Real       | Imag       | Mag        | Expected Real | Expected Imag | Expected Mag | Abs Error" << std::endl;
//     std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
//     bool pass=true;
//     for (size_t i = 0; i < N; ++i)    {
//         double r = (double)input[i].real();
//         double img = (double)input[i].imag();
//         double mag = std::sqrt(r * r + img * img);
//         double exp_r = (double)expected_result[i].real();
//         double exp_img = (double)expected_result[i].imag();
//         double exp_mag = std::sqrt(exp_r * exp_r + exp_img * exp_img);
//         double abs_error = std::sqrt((r - exp_r) * (r - exp_r) + (img - exp_img) * (img - exp_img));
//         if (abs_error > 1e-4) {
//             pass = false;
//         }
//         std::cout << std::setw(5) << i << " | "
//                   << std::fixed << std::setprecision(5)
//                   << std::setw(10) << r << " | "
//                   << std::setw(10) << img << " | "
//                   << std::setw(10) << mag << " | "
//                   << std::setw(13) << exp_r << " | "
//                   << std::setw(13) << exp_img << " | "
//                   << std::setw(13) << exp_mag << " | "
//                   << std::setw(9) << abs_error << std::endl;
//     }
//     if (pass) {
//         std::cout << "✨ Test passed!" << std::endl;
//     } else {
//         std::cout << "❌ Test failed!" << std::endl;
//     }
// }

int main()
{
    test_known_complex_sequence<Q23,Q23Complex>();
    test_known_complex_sequence<double, std::complex<double>>();
    // test_ifft_known_complex_sequence<Q23Complex>();
    // test_ifft_known_complex_sequence<std::complex<double>>();
    return 0;
}