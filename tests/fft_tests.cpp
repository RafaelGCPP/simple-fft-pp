#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "fft_dif.h"
#include "fixed_point.h"

void test_impulse_response()
{
    const size_t N = 16;
    std::cout << "--- Testing FFT Forward DIF (Impulse Response N=" << N << ") ---" << std::endl;

    // 1. Setup Signal: Impulse at index 0
    // We use Q23Complex for the signal buffer
    std::vector<Q23Complex> buffer(N, Q23Complex(0.0, 0.0));
    buffer[0] = Q23Complex(1.0, 0.0);

    // 2. Setup FFT
    using MyTwidGen = TwiddleGenerator<Q31Complex, N>;
    using MyFFT = ForwardFFT_DIF<Q23Complex, MyTwidGen>;

    // 3. Execute
    MyFFT::process(buffer.data());

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

    if (all_ok)
    {
        std::cout << "\n✨ FFT Impulse Test Passed!" << std::endl;
    }
    else
    {
        std::cout << "\n⚠️ FFT Impulse Test Failed (Check scaling/butterfly)." << std::endl;
    }
}

// Helper to calculate magnitude
double get_mag(Q23Complex c)
{
    double r = (double)c.real();
    double i = (double)c.imag();
    return std::sqrt(r * r + i * i);
}

void test_sine_wave()
{
    const size_t N = 16;
    const int target_bin = 1; // 1 cycle per window
    std::cout << "\n--- Testing FFT Forward DIF (Sine Wave k=" << target_bin << ", N=" << N << ") ---" << std::endl;

    std::vector<Q23Complex> buffer(N);

    // 1. Generate Sine Wave: x[n] = 0.5 * sin(2 * pi * target_bin * n / N)
    // We use 0.5 amplitude to keep headroom safe
    for (size_t n = 0; n < N; ++n)
    {
        double angle = 2.0 * M_PI * target_bin * n / N;
        buffer[n] = Q23Complex(std::sin(angle), 0.0);
    }

    // 2. Setup FFT
    using MyTwidGen = TwiddleGenerator<Q31Complex, N>;
    using MyFFT = ForwardFFT_DIF<Q23Complex, MyTwidGen>;

    // 3. Execute
    MyFFT::process(buffer.data());

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
}

int main()
{
    test_impulse_response();
    test_sine_wave();
    return 0;
}