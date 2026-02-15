#pragma once

#include <complex>
#include "fft_core.h"

namespace sfft {  

using Complex = std::complex<float>;

// Shortcuts for Complex FFTs with different sizes
using ComplexFFT1024 = FFT<Complex, Complex,1024>;
using ComplexFFT512 = FFT<Complex, Complex,512>;
using ComplexFFT256 = FFT<Complex, Complex,256>;
using ComplexFFT128 = FFT<Complex, Complex,128>;
using ComplexFFT64 = FFT<Complex, Complex,64>;
using ComplexFFT32 = FFT<Complex, Complex,32>;
using ComplexFFT16 = FFT<Complex, Complex,16>;


} // namespace sfft

