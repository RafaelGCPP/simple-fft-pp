#pragma once

#include <complex>
#include "fft_core.h"
#include "rfft.h"

namespace sfft
{

    using Complex = std::complex<float>;
    using ComplexDouble = std::complex<double>;

    // Shortcuts for Complex FFTs with different sizes
    using ComplexFFT1024 = FFT<Complex, Complex, 1024>;
    using ComplexFFT512 = FFT<Complex, Complex, 512>;
    using ComplexFFT256 = FFT<Complex, Complex, 256>;
    using ComplexFFT128 = FFT<Complex, Complex, 128>;
    using ComplexFFT64 = FFT<Complex, Complex, 64>;
    using ComplexFFT32 = FFT<Complex, Complex, 32>;
    using ComplexFFT16 = FFT<Complex, Complex, 16>;

    // Shortcuts for Complex RFFTs with different sizes
    using ComplexRFFT1024 = RFFT<Complex, Complex, Complex, 1024>;
    using ComplexRFFT512 = RFFT<Complex, Complex, Complex, 512>;
    using ComplexRFFT256 = RFFT<Complex, Complex, Complex, 256>;
    using ComplexRFFT128 = RFFT<Complex, Complex, Complex, 128>;
    using ComplexRFFT64 = RFFT<Complex, Complex, Complex, 64>;
    using ComplexRFFT32 = RFFT<Complex, Complex, Complex, 32>;
    using ComplexRFFT16 = RFFT<Complex, Complex, Complex, 16>;

    // Shortcuts for ComplexDouble FFTs with different sizes
    using ComplexDoubleFFT1024 = FFT<ComplexDouble, ComplexDouble, 1024>;
    using ComplexDoubleFFT512 = FFT<ComplexDouble, ComplexDouble, 512>;
    using ComplexDoubleFFT256 = FFT<ComplexDouble, ComplexDouble, 256>;
    using ComplexDoubleFFT128 = FFT<ComplexDouble, ComplexDouble, 128>;
    using ComplexDoubleFFT64 = FFT<ComplexDouble, ComplexDouble, 64>;
    using ComplexDoubleFFT32 = FFT<ComplexDouble, ComplexDouble, 32>;
    using ComplexDoubleFFT16 = FFT<ComplexDouble, ComplexDouble, 16>;

    // Shortcuts for ComplexDouble RFFTs with different sizes
    using ComplexDoubleRFFT1024 = RFFT<ComplexDouble, ComplexDouble, ComplexDouble, 1024>;
    using ComplexDoubleRFFT512 = RFFT<ComplexDouble, ComplexDouble, ComplexDouble, 512>;
    using ComplexDoubleRFFT256 = RFFT<ComplexDouble, ComplexDouble, ComplexDouble, 256>;
    using ComplexDoubleRFFT128 = RFFT<ComplexDouble, ComplexDouble, ComplexDouble, 128>;
    using ComplexDoubleRFFT64 = RFFT<ComplexDouble, ComplexDouble, ComplexDouble, 64>;
    using ComplexDoubleRFFT32 = RFFT<ComplexDouble, ComplexDouble, ComplexDouble, 32>;
    using ComplexDoubleRFFT16 = RFFT<ComplexDouble, ComplexDouble, ComplexDouble, 16>;

} // namespace sfft
