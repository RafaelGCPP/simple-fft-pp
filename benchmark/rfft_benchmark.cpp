#include <iostream>
#include <complex>
#include <utility>
#include "simple_fft.h"
#include "rfft.h"

template <typename T>
T generate_sample()
{
    return T((double)(random() - 0x40000000) / 0x80000000l);
}

template <typename T, typename CplxT, typename TwidT, int N>
void rfft_benchmark_size(T *data)
{
    std::cout << N << "-points real-valued FFT: ";

    using FullTwidGen = sfft::TwiddleGenerator<TwidT, N>;

    auto fft = sfft::RFFT<T, CplxT, FullTwidGen>();
    auto ifft = sfft::IRFFT<T, CplxT, FullTwidGen>();

    unsigned int start = clock();
    for (int i = 0; i < 5000; i++)
    {
        fft.process(data);
        ifft.process(data);
    }
    unsigned int elapsed = clock() - start;
    std::cout << elapsed / 10000.0 << "us per transform" << std::endl;
}

template <typename T, typename CplxT, typename TwidT, int... Ns>
void rfft_benchmark_sizes(T *data, std::integer_sequence<int, Ns...>)
{
    (rfft_benchmark_size<T, CplxT, TwidT, Ns>(data), ...);
}

template <typename T, typename CplxT, typename TwidT>
void rfft_benchmark()
{
    std::cout << "-=-=-=-=-=-= Real FFT benchmark - "<< typeid(T).name() <<" =-=-=-=-=-=-" << std::endl;

    T data[2048];

    std::cout << "Generating sample data." << std::endl;
    for (int i = 0; i < 2048; i++)
    {
        data[i] = generate_sample<T>();
    }

    rfft_benchmark_sizes<T, CplxT, TwidT>(data,
        std::integer_sequence<int, 8, 16, 32, 64, 128, 256, 512, 1024, 2048>{});
}

int main()
{
    // rfft_benchmark<float, std::complex<float>, std::complex<float>>();
    rfft_benchmark<double, std::complex<double>, std::complex<double>>();
    // rfft_benchmark<sfft::Q23, sfft::Q23Complex, sfft::Q31Complex>();
    return 0;
}