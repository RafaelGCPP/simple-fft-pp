#pragma once

#include "fixed_point.h"
#include "fft_core.h"
#include "rfft.h"

namespace sfft
{

    // Shortcuts for Fixed point Q23 FFTs with different sizes
    using FixedFFT1024Q23 = FFT<Q23Complex, Q31Complex, 1024>;
    using FixedFFT512Q23 = FFT<Q23Complex, Q31Complex, 512>;
    using FixedFFT256Q23 = FFT<Q23Complex, Q31Complex, 256>;
    using FixedFFT128Q23 = FFT<Q23Complex, Q31Complex, 128>;
    using FixedFFT64Q23 = FFT<Q23Complex, Q31Complex, 64>;
    using FixedFFT32Q23 = FFT<Q23Complex, Q31Complex, 32>;
    using FixedFFT16Q23 = FFT<Q23Complex, Q31Complex, 16>;

    // Shortcuts for Fixed point Q15 FFTs with different sizes
    using FixedFFT1024Q15 = FFT<Q15Complex, Q31Complex, 1024>;
    using FixedFFT512Q15 = FFT<Q15Complex, Q31Complex, 512>;
    using FixedFFT256Q15 = FFT<Q15Complex, Q31Complex, 256>;
    using FixedFFT128Q15 = FFT<Q15Complex, Q31Complex, 128>;
    using FixedFFT64Q15 = FFT<Q15Complex, Q31Complex, 64>;
    using FixedFFT32Q15 = FFT<Q15Complex, Q31Complex, 32>;
    using FixedFFT16Q15 = FFT<Q15Complex, Q31Complex, 16>;

    // Shortcuts for Fixed point Q23 RFFTs with different sizes
    using FixedRFFT1024Q23 = RFFT<Q23Complex, Q23Complex, Q31Complex, 1024>;
    using FixedRFFT512Q23 = RFFT<Q23Complex, Q23Complex, Q31Complex, 512>;
    using FixedRFFT256Q23 = RFFT<Q23Complex, Q23Complex, Q31Complex, 256>;
    using FixedRFFT128Q23 = RFFT<Q23Complex, Q23Complex, Q31Complex, 128>;
    using FixedRFFT64Q23 = RFFT<Q23Complex, Q23Complex, Q31Complex, 64>;
    using FixedRFFT32Q23 = RFFT<Q23Complex, Q23Complex, Q31Complex, 32>;
    using FixedRFFT16Q23 = RFFT<Q23Complex, Q23Complex, Q31Complex, 16>;

    // Shortcuts for Fixed point Q15 RFFTs with different sizes
    using FixedRFFT1024Q15 = RFFT<Q15Complex, Q15Complex, Q31Complex, 1024>;
    using FixedRFFT512Q15 = RFFT<Q15Complex, Q15Complex, Q31Complex, 512>;
    using FixedRFFT256Q15 = RFFT<Q15Complex, Q15Complex, Q31Complex, 256>;
    using FixedRFFT128Q15 = RFFT<Q15Complex, Q15Complex, Q31Complex, 128>;
    using FixedRFFT64Q15 = RFFT<Q15Complex, Q15Complex, Q31Complex, 64>;
    using FixedRFFT32Q15 = RFFT<Q15Complex, Q15Complex, Q31Complex, 32>;
    using FixedRFFT16Q15 = RFFT<Q15Complex, Q15Complex, Q31Complex, 16>;

} // namespace sfft