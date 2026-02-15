#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include <numbers> // for std::numbers::pi
#include <simple_fft.h>

using namespace sfft;

int test_impulse_response()
{
    const size_t N = 16;
    std::cout << "--- Testing FFT Forward DIF (Impulse Response N=" << N << ") ---" << std::endl;

    // 1. Setup Signal: Impulse at index 0
    // We use Q23Complex for the signal buffer
    std::vector<Q23Complex> buffer(N, Q23Complex(0.0, 0.0));
    buffer[0] = Q23Complex(1.0, 0.0);

    // 2. Setup FFT
    auto fft = FixedFFT16Q23(); // FFT<Q23Complex, Q31Complex, N>();

    // 3. Execute
    fft.process(buffer.data());

    // 4. Verify
    // Since input was (1,0,0...), output bins should all have magnitude 1.0
    // (Note: Bit-reversal doesn't matter for magnitude in an impulse test)
    bool all_ok = true;
    std::cout << "Index | Real     | Imag     | Mag" << std::endl;
    std::cout << "-----------------------------------" << std::endl;

    for (size_t i = 0; i < N; ++i)
    {
        double r = (double)buffer[i].real();
        double img = (double)buffer[i].imag();
        double mag = std::sqrt(r * r + img * img);

        std::cout << std::setw(5) << i << " | "
                  << std::fixed << std::setprecision(4)
                  << std::setw(8) << r << " | "
                  << std::setw(8) << img << " | "
                  << std::setw(8) << mag;

        if (std::abs(mag - 1.0) > 0.01)
        {
            std::cout << " ❌";
            all_ok = false;
        }
        else
        {
            std::cout << " ✅";
        }
        std::cout << std::endl;
    }

    if (!all_ok)

    {
        std::cout << "\n⚠️ FFT Impulse Test Failed (Check scaling/butterfly)." << std::endl;
        return -1;
    }

    std::cout << "\n✨ FFT Impulse Test Passed!" << std::endl;
    return 0;
}

// Helper to calculate magnitude
double get_mag(Q23Complex c)
{
    double r = (double)c.real();
    double i = (double)c.imag();
    return std::sqrt(r * r + i * i);
}

int test_sine_wave()
{
    const size_t N = 16;
    const int target_bin = 1; // 1 cycle per window
    std::cout << "\n--- Testing FFT Forward DIF (Sine Wave k=" << target_bin << ", N=" << N << ") ---" << std::endl;

    std::vector<Q23Complex> buffer(N);

    // 1. Generate Sine Wave: x[n] = 0.5 * sin(2 * pi * target_bin * n / N)
    // We use 0.5 amplitude to keep headroom safe
    for (size_t n = 0; n < N; ++n)
    {
        double angle = 2.0 * std::numbers::pi * target_bin * n / N;
        buffer[n] = Q23Complex(std::sin(angle), 0.0);
    }

    // 2. Setup FFT
    auto fft = FixedFFT16Q23(); // FFT<Q23Complex, Q31Complex, N>();

    // 3. Execute
    fft.process(buffer.data());

    // 4. Verify results
    // For a sine wave, we expect peaks at bin 1 and bin 15 (N-1)
    // Note: Since this is DIF, the output is BIT-REVERSED.
    // Index 1 bit-reversed (4 bits) is 1000 binary = 8
    // Index 15 bit-reversed (4 bits) is 1111 binary = 15

    std::cout << std::endl;
    std::cout << "Using raw data" << std::endl;

    std::cout << "Raw Idx | Real     | Imag     | Mag" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    for (size_t i = 0; i < N; ++i)
    {
        double mag = get_mag(buffer[i]);

        std::cout << std::setw(7) << i << " | "
                  << std::fixed << std::setprecision(4)
                  << std::setw(8) << (double)buffer[i].real() << " | "
                  << std::setw(8) << (double)buffer[i].imag() << " | "
                  << std::setw(8) << mag;

        if (mag > 1.0)
            std::cout << " 🔥 PEAK";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Using view decorator" << std::endl;

    std::cout << "Idx     | Real     | Imag     | Mag" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    BitReversedView<Q23Complex, N> view(buffer.data());
    for (size_t i = 0; i < N; ++i)
    {
        double mag = get_mag(view[i]);

        std::cout << std::setw(7) << i << " | "
                  << std::fixed << std::setprecision(4)
                  << std::setw(8) << (double)view[i].real() << " | "
                  << std::setw(8) << (double)view[i].imag() << " | "
                  << std::setw(8) << mag;

        if (mag > 1.0)
            std::cout << " 🔥 PEAK";
        std::cout << std::endl;
    }

    view.commit(); // Now the buffer is in normal order
    std::cout << std::endl;
    std::cout << "Using committed data" << std::endl;

    std::cout << "Raw Idx | Real     | Imag     | Mag" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    for (size_t i = 0; i < N; ++i)
    {
        double mag = get_mag(buffer[i]);

        std::cout << std::setw(7) << i << " | "
                  << std::fixed << std::setprecision(4)
                  << std::setw(8) << (double)buffer[i].real() << " | "
                  << std::setw(8) << (double)buffer[i].imag() << " | "
                  << std::setw(8) << mag;

        if (mag > 1.0)
            std::cout << " 🔥 PEAK";
        std::cout << std::endl;
    }
    return 0;
}

template <typename T, typename U>
int test_known_complex_sequence()
{
    const size_t N = 8;
    std::cout << "\n--- Testing FFT Forward DIF (Known Sequence N=" << N << ") ---" << std::endl;
    T input[N] = {
        T(1.0, 2.0),
        T(3.0, 4.0),
        T(5.0, 6.0),
        T(7.0, 8.0),
        T(-8.0, -7.0),
        T(-6.0, -5.0),
        T(-4.0, -3.0),
        T(-2.0, -1.0)};

    T expected_fft[N] = {
        T(-4.0, 4.0),
        T(30.72792, -12.72792),
        T(-16.0, 0.0),
        T(12.72792, 5.27208),
        T(-8.0, -8.0),
        T(5.27208, 12.72792),
        T(0.0, -16.0),
        T(-12.72792, 30.72792)};

    auto fft = FFT<T, U, N>();
    auto view = BitReversedView<T, N>(input);
    fft.process(input);
    view.commit(); // Ensure data is in normal order for comparison
    std::cout << "Index | Real       | Imag       | Mag        | Expected Real | Expected Imag | Expected Mag | Abs Error" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    bool pass = true;
    for (size_t i = 0; i < N; ++i)
    {
        double r = (double)input[i].real();
        double img = (double)input[i].imag();
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
    if (!pass)
    {
        std::cout << "❌ Test failed!" << std::endl;
        return -1;
    }

    std::cout << "✨ Test passed!" << std::endl;
    return 0;
}

template <typename T, typename U>
int test_ifft_known_complex_sequence()
{
    const size_t N = 8;
    std::cout << "\n--- Testing IFFT (Known Sequence N=" << N << ") ---" << std::endl;
    T expected_result[N] = {
        T(1.0, 2.0),
        T(3.0, 4.0),
        T(5.0, 6.0),
        T(7.0, 8.0),
        T(-8.0, -7.0),
        T(-6.0, -5.0),
        T(-4.0, -3.0),
        T(-2.0, -1.0)};

    T input[N] = {
        T(-4.0, 4.0),
        T(30.72792, -12.72792),
        T(-16.0, 0.0),
        T(12.72792, 5.27208),
        T(-8.0, -8.0),
        T(5.27208, 12.72792),
        T(0.0, -16.0),
        T(-12.72792, 30.72792)};

    auto fft = FFT<T, U, N>();
    auto view = BitReversedView<T, N>(input);
    view.commit(); // Ensure data is bitreversed for IFFT

    fft.inverse(input);
    std::cout << "Index | Real       | Imag       | Mag        | Expected Real | Expected Imag | Expected Mag | Abs Error" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
    bool pass = true;
    for (size_t i = 0; i < N; ++i)
    {
        double r = (double)input[i].real();
        double img = (double)input[i].imag();
        double mag = std::sqrt(r * r + img * img);
        double exp_r = (double)expected_result[i].real();
        double exp_img = (double)expected_result[i].imag();
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
    if (!pass)
    {
        std::cout << "❌ Test failed!" << std::endl;
        return -1;
    }
    std::cout << "✨ Test passed!" << std::endl;
    return 0;
}

int main()
{
    int retval=0;
    retval |= test_impulse_response();
    retval |= test_sine_wave();
    retval |= test_known_complex_sequence<Q23Complex, Q31Complex>();
    retval |= test_known_complex_sequence<std::complex<double>, std::complex<double>>();
    retval |= test_ifft_known_complex_sequence<Q23Complex, Q31Complex>();
    retval |= test_ifft_known_complex_sequence<std::complex<double>, std::complex<double>>();
    return retval;
}