#include <iostream>
#include <complex>
#include <utility>
#include "simple_fft.h"

template <typename T>
T generate_complex_sample()
{
    return T(((double)(random() - 0x40000000) / 0x80000000l), ((double)(random() - 0x40000000) / 0x80000000l));
}

template<typename CplxT, typename TwidT>
void fft_benchmark()
{
    std::cout << "-=-=-=-=-=-= Complex FFT benchmark =-=-=-=-=-=-" << std::endl;

    CplxT data[2048];

    std::cout << "Generating sample data." << std::endl;
    for (int i = 0; i < 2048; i++)
    {
        data[i] = generate_complex_sample<CplxT>();
    }

    auto fft = sfft::FFT<CplxT, TwidT, 2048>();

    unsigned int start = clock();
    for (int i = 0; i < 5000; i++)
    {
        fft.process(data);
        fft.inverse(data);
    }
    unsigned int elapsed = clock() - start;
    std::cout << elapsed / 10000.0 << "us per operation (forward or inverse)" << std::endl;
}

int main()
{
    fft_benchmark<std::complex<float>, std::complex<float>>();
    fft_benchmark<std::complex<double>, std::complex<double>>();
    fft_benchmark<sfft::Q23Complex, sfft::Q31Complex>();
    return 0;
}