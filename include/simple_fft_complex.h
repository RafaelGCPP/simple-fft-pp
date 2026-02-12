#pragma once

#include <complex>
#include "fft_core.h"

using Complex=std::complex<double>;
using TwiddleGenerator1024 = TwiddleGenerator<Complex, 1024>;
using TwiddleGenerator512 = TwiddleGenerator<Complex, 512>;
using TwiddleGenerator256 = TwiddleGenerator<Complex, 256>;
using TwiddleGenerator128 = TwiddleGenerator<Complex, 128>;
using TwiddleGenerator64 = TwiddleGenerator<Complex, 64>;
using TwiddleGenerator32 = TwiddleGenerator<Complex, 32>;
using TwiddleGenerator16 = TwiddleGenerator<Complex, 16>;    

// Shortcuts for Complex FFTs with different sizes
using ComplexFFT1024 = FFT<Complex, TwiddleGenerator1024>;
using ComplexFFT512 = FFT<Complex, TwiddleGenerator512>;
using ComplexFFT256 = FFT<Complex, TwiddleGenerator256>;
using ComplexFFT128 = FFT<Complex, TwiddleGenerator128>;
using ComplexFFT64 = FFT<Complex, TwiddleGenerator64>;
using ComplexFFT32 = FFT<Complex, TwiddleGenerator32>;
using ComplexFFT16 = FFT<Complex, TwiddleGenerator16>;


// Shortcuts for Complex IFFTs with different sizes
using ComplexIFFT1024 = IFFT<Complex, TwiddleGenerator1024>;
using ComplexIFFT512 = IFFT<Complex, TwiddleGenerator512>;
using ComplexIFFT256 = IFFT<Complex, TwiddleGenerator256>;
using ComplexIFFT128 = IFFT<Complex, TwiddleGenerator128>;
using ComplexIFFT64 = IFFT<Complex, TwiddleGenerator64>;
using ComplexIFFT32 = IFFT<Complex, TwiddleGenerator32>;
using ComplexIFFT16 = IFFT<Complex, TwiddleGenerator16>;

