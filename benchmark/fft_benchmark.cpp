#include <iostream>
#include <complex>
#include <utility>
#include <random>
#include <ctime>
#include "simple_fft.h"

std::mt19937 gen(std::random_device{}());
std::uniform_real_distribution<double> dist(-0.5, 0.5);
constexpr double timescale_us = 1000000.0 / CLOCKS_PER_SEC;

template <typename T>
T generate_complex_sample()
{
    return T(dist(gen), dist(gen));
}

template <typename CplxT, typename TwidT, int N>
void fft_benchmark_size(CplxT *data)
{
    std::cout << N << "-points complex FFT: ";

    auto fft = sfft::FFT<CplxT, TwidT, N>();

    unsigned int start = clock();
    for (int i = 0; i < 5000; i++)
    {
        fft.process(data);
        fft.inverse(data);
    }
    unsigned int elapsed = clock() - start;
    std::cout << elapsed * timescale_us / 10000.0 << "us per transform" << std::endl;
}

template <typename CplxT, typename TwidT, int... Ns>
void fft_benchmark_sizes(CplxT *data, std::integer_sequence<int, Ns...>)
{
    (fft_benchmark_size<CplxT, TwidT, Ns>(data), ...);
}

template <typename CplxT, typename TwidT>
void fft_benchmark()
{
    std::cout << "-=-=-=-=-=-= Complex FFT benchmark - " << typeid(CplxT).name() << " =-=-=-=-=-=-" << std::endl;

    CplxT data[2048];

    std::cout << "Generating sample data." << std::endl;
    for (int i = 0; i < 2048; i++)
    {
        data[i] = generate_complex_sample<CplxT>();
    }

    fft_benchmark_sizes<CplxT, TwidT>(data,
        std::integer_sequence<int, 8, 16, 32, 64, 128, 256, 512, 1024, 2048>{});
}

int main()
{
    fft_benchmark<std::complex<float>, std::complex<float>>();
    fft_benchmark<std::complex<double>, std::complex<double>>();
    fft_benchmark<sfft::Q23Complex, sfft::Q31Complex>();
    return 0;
}