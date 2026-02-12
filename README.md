# Simple FFT in C++

A header-only C++20 meta-programming library for computing Fast Fourier Transforms (FFT) with compile-time unrolling and fixed-point arithmetic support. Designed for embedded systems and microcontroller platforms.

## Features

- **Header-only library** — no compilation needed for the FFT logic itself
- **Compile-time unrolled loops** — zero-overhead abstractions via C++20 constexpr
- **Mixed-precision arithmetic** — Q15, Q23, Q31 fixed-point types with cross-type operations
- **DIF/DIT asymmetry** — Forward FFT uses Decimation in Frequency (DIF), Inverse uses Decimation in Time (DIT), eliminating the need for input bit-reversal
- **Bit-reversal decorator** — optional in-place reordering with lazy evaluation
- **Microcontroller-friendly** — avoids features incompatible with RP2350 and similar platforms

## Building

### Prerequisites

- C++20 compiler (GCC 10+, Clang 10+, MSVC 2019+)
- CMake 3.16+

### Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```

### Run Tests

```bash
./fft_tests
./roundtrip_tests
./fixed_point_tests
./butterfly_tests
./twiddles_test
```

## Library Usage

### Basic Example

Include the main header:

```cpp
#include "simple_fft.h"

// For double precision
using TwidGen = TwiddleGenerator<std::complex<double>, 64>;
using FFT = ForwardFFT_DIF<std::complex<double>, TwidGen>;

int main()
{
    std::vector<std::complex<double>> data(64);
    // ... populate data ...
    FFT::process(data.data());
    return 0;
}
```

### Fixed-Point Example

```cpp
#include "simple_fft.h"

// Use Q15 (15 fractional bits) for data, Q31 for twiddles (higher precision)
using TwidGen = TwiddleGenerator<Q31Complex, 64>;
using FFT = ForwardFFT_DIF<Q15Complex, TwidGen>;

int main()
{
    std::vector<Q15Complex> data(64);
    // ... populate data ...
    FFT::process(data.data());
    return 0;
}
```

### Round-trip (FFT → IFFT)

```cpp
#include "simple_fft.h"

int main()
{
    using TwidGen = TwiddleGenerator<std::complex<double>, 64>;
    using FwdFFT = ForwardFFT_DIF<std::complex<double>, TwidGen>;
    using InvFFT = InverseFFT_DIT<std::complex<double>, TwidGen>;

    std::vector<std::complex<double>> signal(64);
    // ... populate signal ...

    FwdFFT::process(signal.data());  // Forward FFT
    InvFFT::process(signal.data());  // Inverse FFT
    // signal is now restored (approximately, modulo numerical errors)
    
    return 0;
}
```

## Real-World Scenarios

### Audio Signal Processing on Microcontrollers
Analyze audio chunks for equalization or visualization.
```cpp
// 128-point FFT for a small audio buffer
// Input: Q15 audio samples (e.g. from a microphone ADC)
using AudioFFT = ForwardFFT_DIF<Q15Complex, TwiddleGenerator<Q31Complex, 128>>;

void process_audio_chunk(Q15Complex* buffer) {
    AudioFFT::process(buffer);
    // buffer now contains frequency bins in bit-reversed order.
    // Useful for spectral display or detecting dominant frequencies.
}
```

### Frequency Domain Filtering
Efficiently apply filters in the frequency domain. Since the Forward FFT produces bit-reversed output, we can avoid runtime reordering by pre-calculating the filter mask in bit-reversed order.
```cpp
// 256-point FFT/IFFT for filtering
using FwdFFT = ForwardFFT_DIF<std::complex<double>, TwiddleGenerator<std::complex<double>, 256>>;
using InvFFT = InverseFFT_DIT<std::complex<double>, TwiddleGenerator<std::complex<double>, 256>>;

// Pre-calculated filter coefficients stored in BIT-REVERSED order
extern const std::vector<std::complex<double>> bit_reversed_filter_mask;

void apply_filter(std::vector<std::complex<double>>& signal) {
    // 1. Transform to frequency domain (Output is bit-reversed)
    FwdFFT::process(signal.data());

    // 2. Apply filter mask (Element-wise multiplication)
    // No bit-reversal needed because the mask matches the FFT output order
    for(size_t i = 0; i < 256; ++i) {
        signal[i] *= bit_reversed_filter_mask[i];
    }

    // 3. Transform back to time domain (Input bit-reversed -> Output natural)
    InvFFT::process(signal.data());
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

std::vector<std::complex<double>> data(64);
// ... populate data ...

// Wrap data in a decorator for lazy bit-reversal access
BitReversedView<std::complex<double>, 64> view(data.data());

// Access element at logical index 5 (physically bit-reversed)
auto element = view[5];  // Returns data at bit-reversed index

// Commit the bit-reversal: physically reorder data in-place
view.commit();  // Now data is physically bit-reversed

// Subsequent access is O(1) without reversals
auto element2 = data[5];  // Direct access to reordered data
```

**Use Cases**:

- **Sparse Access / Lazy Evaluation**:
  Access specific frequency bins by their natural index without reordering the entire array. Useful when you only need to check a few dominant frequencies.

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

*Note: The underlying storage is always `int32_t`. Adjust `B` to trade off between integer range (bits = 31 - B) and fractional precision.*

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
simple_fft.h
├── fft_core.h          — Butterfly operations, bit-reversal utilities
├── twiddles.h          — Twiddle factor generation (exp(-2πj*k/N))
├── fixed_point.h       — Q15/Q23/Q31 types with mixed-precision arithmetic
├── fft_dif.h           — Forward FFT (Decimation in Frequency)
└── ifft_dit.h          — Inverse FFT (Decimation in Time)
```

## Compiler Support

- **Recommended**: GCC 10+, Clang 12+
- **C++ Standard**: C++20 (C++17 minimum with some features disabled)
- **Platforms**: x86_64, ARM (incl. RP2350), MIPS, others with C++20 support

## Testing

The library includes comprehensive test suites:

- `fft_tests.cpp` — Impulse and frequency-domain FFT verification
- `roundtrip_tests.cpp` — Signal reconstruction (FFT → IFFT)
- `fixed_point_tests.cpp` — Arithmetic correctness for fixed-point types
- `butterfly_tests.cpp` — Individual butterfly operation verification
- `twiddles_test.cpp` — Twiddle factor accuracy

Run with:
```bash
cmake .. && make && make test
```

## License

[Add your license here]

## Contributing

Contributions are welcome! Please ensure:
- All code is in English (comments, documentation, commit messages)
- Code targets C++20 where possible (C++17 for microcontroller compatibility)
- All tests pass before submitting PRs
