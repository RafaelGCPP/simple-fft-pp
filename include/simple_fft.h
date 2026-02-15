#pragma once

#include "fixed_point.h"
#include "fft_core.h"

namespace sfft {

// Shortcuts for Fixed point Q23 FFTs with different sizes
using FixedFFT1024Q23 = FFT<Q23Complex, Q31Complex,1024>;
using FixedFFT512Q23 = FFT<Q23Complex, Q31Complex,512>;
using FixedFFT256Q23 = FFT<Q23Complex, Q31Complex,256>;
using FixedFFT128Q23 = FFT<Q23Complex, Q31Complex,128>;
using FixedFFT64Q23 = FFT<Q23Complex, Q31Complex,64>;
using FixedFFT32Q23 = FFT<Q23Complex, Q31Complex,32>;
using FixedFFT16Q23 = FFT<Q23Complex, Q31Complex,16>;

//Shortcuts for Fixed point Q15 FFTs with different sizes
using FixedFFT1024Q15 = FFT<Q15Complex, Q31Complex,1024>;
using FixedFFT512Q15 = FFT<Q15Complex, Q31Complex,512>;
using FixedFFT256Q15 = FFT<Q15Complex, Q31Complex,256>;
using FixedFFT128Q15 = FFT<Q15Complex, Q31Complex,128>;
using FixedFFT64Q15 = FFT<Q15Complex, Q31Complex,64>;
using FixedFFT32Q15 = FFT<Q15Complex, Q31Complex,32>;
using FixedFFT16Q15 = FFT<Q15Complex, Q31Complex,16>;

} // namespace sfft