# Simple FFT++

A **header-only C++20 library** that generates fully unrolled, statically dispatched FFT routines at compile time. It targets **embedded and microcontroller platforms** (RP2350, ARM Cortex-M, etc.) where dynamic memory allocation and runtime planning are unacceptable.

The library provides:

- **Complex FFT / IFFT** — arbitrary power-of-2 sizes, both floating-point and fixed-point
- **Real FFT / Inverse RFFT** — computes the spectrum of real-valued signals at half the cost, with an in-place spectral manipulation API
- **Fixed-point arithmetic** — Q15, Q23, Q31 types with mixed-precision operations, no floats required at runtime

> **This is not a general-purpose FFT library.** It is designed for scenarios where the FFT size is known at compile time and the priority is deterministic, allocation-free execution on resource-constrained hardware.

## Features

- **Header-only** — drop the `include/` folder into your project, no build step for the library itself
- **Compile-time unrolling** — all $\log_2 N$ stages are expanded by the compiler into a flat call sequence
- **Zero dynamic allocation** — everything operates in-place on caller-provided buffers
- **Mixed-precision fixed-point** — signal data in Q15/Q23 with twiddle factors in Q31 for maximum accuracy
- **DIF/DIT asymmetry** — Forward FFT (DIF) accepts natural-order input, Inverse FFT (DIT) produces natural-order output — no explicit bit-reversal pass needed for round-trips
- **Real FFT with lazy unpacking** — `RFFT_View` provides O(1) access to any spectral bin without materializing the full spectrum
- **In-place spectral transforms** — modify frequency bins through a lambda and reconstruct, all without auxiliary buffers
- **Microcontroller-friendly** — avoids exceptions, RTTI, heap allocation, and virtual dispatch

## Quick Start

### Complex FFT (Floating-Point)

```cpp
#include "simple_fft_complex.h"   // convenience aliases for std::complex<float>
using namespace sfft;

int main() {
    Complex data[64];
    // ... populate data ...

    auto fft = ComplexFFT64();
    fft.process(data);    // Forward FFT (output in bit-reversed order)
    fft.inverse(data);    // Inverse FFT (output in natural order)
    // data is restored
}
```

### Complex FFT (Fixed-Point)

```cpp
#include "simple_fft.h"   // convenience aliases for fixed-point types
using namespace sfft;

int main() {
    Q23Complex data[128];
    // ... populate data ...

    auto fft = FixedFFT128Q23();
    fft.process(data);    // Forward FFT
    fft.inverse(data);    // Inverse FFT
    // data is restored
}
```

### Round-trip (FFT → IFFT)

```cpp
#include "fft_core.h"
using namespace sfft;

int main()
{
    using ComplexType = std::complex<double>;
    auto fft = FFT<ComplexType, ComplexType, 64>();

    ComplexType signal[64];
    // ... populate signal ...

    fft.process(signal);  // Forward FFT (output in bit-reversed order)
    fft.inverse(signal);  // Inverse FFT (output in natural order)
    // signal is now restored (approximately, modulo numerical errors)
    
    return 0;
}
```

### Real FFT with Spectral Filtering

For real-valued signals, `RFFT` computes the spectrum at half the cost. See [RealFFT.md](RealFFT.md) for the full mathematical derivation.

```cpp
#include "rfft.h"
using namespace sfft;

int main() {
    const size_t N = 256;
    using ComplexType = std::complex<double>;
    using RealType = ComplexType::value_type;
    
    RealType buffer[N];
    // ... fill buffer with real-valued samples ...

    auto rfft = RFFT<RealType, ComplexType, ComplexType, N>();

    // Forward: N real samples → N/2+1 unique spectral bins
    auto view = rfft.process(buffer);

    // Access any bin in O(1)
    auto dc = view[0];
    auto bin5 = view[5];

    // In-place low-pass filter: zero out bins above cutoff
    view.transform([](size_t k, ComplexType val) {
        if (k > 20) return ComplexType(0.0, 0.0);
        return val;
    });

    // Inverse: reconstruct filtered signal in-place
    rfft.inverse(buffer);
    // buffer now contains the filtered real signal
}
```

## Real-World Scenarios

### Audio Signal Processing on Microcontrollers
Analyze audio chunks for equalization or visualization.
```cpp
// 128-point FFT for a small audio buffer
// Input: Q15 audio samples (e.g. from a microphone ADC)
using namespace sfft;

void process_audio_chunk(Q15Complex* buffer) {
    auto fft = FixedFFT128Q15();
    fft.process(buffer);
    // buffer now contains frequency bins in bit-reversed order.
    // Useful for spectral display or detecting dominant frequencies.
}
```

### Frequency Domain Filtering (Complex)
Efficiently apply filters in the frequency domain. Since the Forward FFT produces bit-reversed output, we can avoid runtime reordering by pre-calculating the filter mask in bit-reversed order.
```cpp
// 256-point FFT/IFFT for filtering
using namespace sfft;

// Pre-calculated filter coefficients stored in BIT-REVERSED order
extern const std::complex<double> bit_reversed_filter_mask[256];

void apply_filter(std::complex<double>* signal) {
    auto fft = ComplexDoubleFFT256();
    
    // 1. Transform to frequency domain (Output is bit-reversed)
    fft.process(signal);

    // 2. Apply filter mask (Element-wise multiplication)
    // No bit-reversal needed because the mask matches the FFT output order
    for(size_t i = 0; i < 256; ++i) {
        signal[i] *= bit_reversed_filter_mask[i];
    }

    // 3. Transform back to time domain (Input bit-reversed -> Output natural)
    fft.inverse(signal);
}
```

### Frequency Domain Filtering (Real)
For real-valued signals, use `RFFT` with the `transform` API — no need for bit-reversed masks or manual repacking. See [RealFFT.md](RealFFT.md) for details.
```cpp
#include "rfft.h"
using namespace sfft;

void apply_lowpass(double* signal) {
    const size_t N = 256;
    using ComplexType = std::complex<double>;
    using RealType = ComplexType::value_type;
    
    auto rfft = RFFT<RealType, ComplexType, ComplexType, N>();

    auto view = rfft.process(signal);

    view.transform([](size_t k, ComplexType val) {
        if (k > 20) return ComplexType(0.0, 0.0);
        return val;
    });

    rfft.inverse(signal);
}
```

## Design Decisions

### Why Metaprogramming?

In embedded real-time systems, predictability and efficiency are paramount. General-purpose FFT libraries often use dynamic planning or runtime recursion, which incurs overhead.

By providing the FFT parameters (size, type, direction) as **template arguments**, this library allows the compiler to:

1. **Calculate Constants at Compile-Time**: Twiddle factor indices, loop strides, and butterfly groups are known constants.
2. **Unroll Recursion**: The $O(\log N)$ stages of the FFT are expanded into a linear sequence of function calls during compilation.
3. **Inline Operations**: Small butterfly operations are easily inlined, removing function call overhead.
4. **Optimize Register Usage**: With constant loop bounds, the compiler can better allocate registers and apply SIMD instructions (like NEON or AVX) automatically.

This results in a "hard-coded" FFT implementation tailored exactly to your specifications, similar to what you might write by hand in assembly, but with C++ type safety and readability.

### DIF for Forward, DIT for Inverse

This library uses an **asymmetric** FFT strategy:

- **Forward FFT**: Decimation in Frequency (DIF)
  - Operates on the **input in natural order**
  - Produces output in **bit-reversed order**
  
- **Inverse FFT**: Decimation in Time (DIT)
  - Operates on input in **bit-reversed order** (as produced by DIF)
  - Produces output in **natural order**

**Benefit**: No need to bit-reverse the input before forward transform. After a round-trip (FFT → IFFT), the signal is naturally ordered without explicit reordering.

```
Signal (natural)                Output (natural)
    ↓                                ↑
  DIF FFT                         DIT IFFT
    ↓                                ↑
Output (bit-reversed)      Input (bit-reversed)
```

### Bit-Reversal Decorator

For cases where you need to manipulate bit-reversed data or perform explicit reordering, the library provides a **`BitReversedView`** decorator:

```cpp
#include "fft_core.h"
using namespace sfft;

std::vector<std::complex<double>> data(64);
// ... populate data ...

// Wrap data in a decorator for lazy bit-reversal access
BitReversedView<std::complex<double>, 64> view(data.data());

// Access element at logical index 5 (physically bit-reversed)
auto element = view[5];  // Returns data at bit-reversed index

// Apply in-place transformations without committing the bit-reversal
view.transform([](size_t k, std::complex<double> val) {
    if (k > 20) return std::complex<double>(0.0, 0.0);  // Zero out high frequencies
    return val;
});

// Commit the bit-reversal: physically reorder data in-place
view.commit();  // Now data is physically bit-reversed

// Subsequent access is O(1) without reversals
auto element2 = data[5];  // Direct access to reordered data
```

**Use Cases**:

- **Sparse Access / Lazy Evaluation**:
  Access specific frequency bins by their natural index without reordering the entire array. Useful when you only need to check a few dominant frequencies.

- **In-Place Spectral Filtering**:
  Apply transformations using the `transform()` method while the data remains in bit-reversed order, avoiding the overhead of reordering twice.

- **Interfacing with Other Libraries**:
  Bridge the gap between this library's bit-reversed output and other algorithms (or visualization tools) that expect data in natural sequence.

- **Spectral Manipulation**:
  Simplify operations like zero-padding (for interpolation) or applying complex frequency masks where reasoning about index mapping is difficult in bit-reversed order.

## Fixed-Point Types

The library uses a generic `FixedPoint<B>` template where `B` specifies the number of fractional bits. While the library pre-defines three common aliases, you can instantiate any Q-format (e.g., `FixedPoint<20>`) as needed.

**Pre-defined Aliases**:

| Type | Alias For | Range | Resolution |
|------|-----------|-------|------------|
| Q15  | `FixedPoint<15>` | [-65536, 65536)   | ~3.05e-5  |
| Q23  | `FixedPoint<23>` | [-256, 256)       | ~1.19e-7  |
| Q31  | `FixedPoint<31>` | [-1, 1)           | ~4.66e-10 |

*Note: The underlying storage is always `int32_t` (32 bits). The integer part uses (31 - B) bits, and the fractional part uses B bits. The range is approximately [−2^(31−B), 2^(31−B)) and the resolution is 2^−B.*

### Complex Counterparts

For FFT operations, scalar types are wrapped in `FixedComplex<T>`, which mimics the interface of `std::complex`.

| Complex Type | Underlying Components |
|--------------|-----------------------|
| `Q15Complex` | `FixedComplex<Q15>`   |
| `Q23Complex` | `FixedComplex<Q23>`   |
| `Q31Complex` | `FixedComplex<Q31>`   |

These types support standard complex arithmetic (conjugate, magnitude squared, multiplication) and can be mixed with `std::complex` for convenience during setup or testing.

Operations between different formats are automatically handled with appropriate scaling.

## Architecture

```
include/
├── fft_core.h              — FFT engine, DIF/DIT butterflies, BitReversedView, TwiddleGenerator
├── constexpr_math.h         — Compile-time sin/cos/pi for twiddle factor generation
├── fixed_point.h           — FixedPoint<B>, FixedComplex<T>, mixed-precision arithmetic
├── rfft.h                  — RFFT, RFFT_View (real-valued FFT)
├── simple_fft.h            — Convenience aliases for fixed-point FFTs (FixedFFT64Q23, etc.)
└── simple_fft_complex.h    — Convenience aliases for std::complex<float/double> FFTs (ComplexFFT64, etc.)
```

## Compiler Support

- **Recommended**: GCC 10+, Clang 12+
- **C++ Standard**: C++20 (C++17 minimum with some features disabled)
- **Platforms**: x86_64, ARM (incl. RP2350), MIPS, others with C++20 support

## Testing

The library includes comprehensive test suites:

- `fft_tests.cpp` — Impulse and frequency-domain FFT verification
- `roundtrip_tests.cpp` — Signal reconstruction (FFT → IFFT)
- `rfft_tests.cpp` — Real FFT correctness, spectral transform round-trip
- `fixed_point_tests.cpp` — Arithmetic correctness for fixed-point types

Run with:
```bash
mkdir build && cd build
cmake ..
make
./tests/fft_tests
./tests/roundtrip_tests
./tests/rfft_tests
./tests/fixed_point_tests
```

## License

[Add your license here]

## Contributing

Contributions are welcome! Please ensure:
- All code is in English (comments, documentation, commit messages)
- Code targets C++20 where possible (C++17 for microcontroller compatibility)
- All tests pass before submitting PRs
