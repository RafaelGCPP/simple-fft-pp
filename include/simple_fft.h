#pragma once

#include "fixed_point.h"
#include "fft_core.h"

namespace sfft {

using TwiddleGenerator1024 = TwiddleGenerator<Q31Complex, 1024>;
using TwiddleGenerator512 = TwiddleGenerator<Q31Complex, 512>;
using TwiddleGenerator256 = TwiddleGenerator<Q31Complex, 256>;
using TwiddleGenerator128 = TwiddleGenerator<Q31Complex, 128>;
using TwiddleGenerator64 = TwiddleGenerator<Q31Complex, 64>;
using TwiddleGenerator32 = TwiddleGenerator<Q31Complex, 32>;
using TwiddleGenerator16 = TwiddleGenerator<Q31Complex, 16>;    

// Shortcuts for Fixed point Q23 FFTs with different sizes
using FixedFFT1024Q23 = FFT<Q23Complex, TwiddleGenerator1024>;
using FixedFFT512Q23 = FFT<Q23Complex, TwiddleGenerator512>;
using FixedFFT256Q23 = FFT<Q23Complex, TwiddleGenerator256>;
using FixedFFT128Q23 = FFT<Q23Complex, TwiddleGenerator128>;
using FixedFFT64Q23 = FFT<Q23Complex, TwiddleGenerator64>;
using FixedFFT32Q23 = FFT<Q23Complex, TwiddleGenerator32>;
using FixedFFT16Q23 = FFT<Q23Complex, TwiddleGenerator16>;

//Shortcuts for Fixed point Q15 FFTs with different sizes
using FixedFFT1024Q15 = FFT<Q15Complex, TwiddleGenerator1024>;
using FixedFFT512Q15 = FFT<Q15Complex, TwiddleGenerator512>;
using FixedFFT256Q15 = FFT<Q15Complex, TwiddleGenerator256>;
using FixedFFT128Q15 = FFT<Q15Complex, TwiddleGenerator128>;
using FixedFFT64Q15 = FFT<Q15Complex, TwiddleGenerator64>;
using FixedFFT32Q15 = FFT<Q15Complex, TwiddleGenerator32>;
using FixedFFT16Q15 = FFT<Q15Complex, TwiddleGenerator16>;

// Shortcuts for Fixed point Q23 IFFTs with different sizes
using FixedIFFT1024Q23 = IFFT<Q23Complex, TwiddleGenerator1024>;
using FixedIFFT512Q23 = IFFT<Q23Complex, TwiddleGenerator512>;
using FixedIFFT256Q23 = IFFT<Q23Complex, TwiddleGenerator256>;
using FixedIFFT128Q23 = IFFT<Q23Complex, TwiddleGenerator128>;
using FixedIFFT64Q23 = IFFT<Q23Complex, TwiddleGenerator64>;
using FixedIFFT32Q23 = IFFT<Q23Complex, TwiddleGenerator32>;
using FixedIFFT16Q23 = IFFT<Q23Complex, TwiddleGenerator16>;

//Shortcuts for Fixed point Q15 IFFTs with different sizes
using FixedIFFT1024Q15 = IFFT<Q15Complex, TwiddleGenerator1024>;
using FixedIFFT512Q15 = IFFT<Q15Complex, TwiddleGenerator512>;
using FixedIFFT256Q15 = IFFT<Q15Complex, TwiddleGenerator256>;
using FixedIFFT128Q15 = IFFT<Q15Complex, TwiddleGenerator128>;
using FixedIFFT64Q15 = IFFT<Q15Complex, TwiddleGenerator64>;
using FixedIFFT32Q15 = IFFT<Q15Complex, TwiddleGenerator32>;
using FixedIFFT16Q15 = IFFT<Q15Complex, TwiddleGenerator16>;

} // namespace sfft