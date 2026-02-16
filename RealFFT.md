# Real FFT (RFFT) — Theory and Implementation

This document explains the mathematical foundation behind the Real FFT and how it is implemented in this library through the `RFFT` and `RFFT_View` classes.

## 1. Motivation

A standard complex FFT of size $N$ operates on $N$ complex inputs and produces $N$ complex outputs. When the input signal is purely real (as is typical in audio, sensor, and control applications), half of the computation is redundant due to the **Hermitian symmetry** of the spectrum:

$$X[k] = \overline{X[N-k]}$$

The Real FFT exploits this by packing $N$ real samples into $N/2$ complex values, computing an $N/2$-point complex FFT, and then recovering the full $N$-point spectrum through an algebraic post-processing step. This halves both memory and computation.

## 2. Mathematical Foundation

### 2.1 Packing Real Data into Complex

Given $N$ real samples $x[0], x[1], \ldots, x[N-1]$, we form $N/2$ complex values:

$$z[n] = x[2n] + j \cdot x[2n+1], \quad n = 0, 1, \ldots, \frac{N}{2} - 1$$

This is equivalent to viewing the real array as interleaved complex pairs — a zero-cost operation via `InterleavedComplexView`.

### 2.2 Relationship Between $Z[k]$ and $X[k]$

Let $Z[k]$ denote the $N/2$-point DFT of the packed sequence $z[n]$, and $X[k]$ the $N$-point DFT of the original real sequence $x[n]$.

Define the **conjugate-symmetric** and **conjugate-antisymmetric** parts of $Z$:

$$A_k = \frac{Z[k] + \overline{Z[\tfrac{N}{2}-k]}}{2}$$

$$B_k = \frac{Z[k] - \overline{Z[\tfrac{N}{2}-k]}}{2}$$

Then the $N$-point DFT bins are recovered as:

$$\boxed{X[k] = A_k - j \cdot B_k \cdot W_N^k}$$

where $W_N^k = e^{-j2\pi k/N}$ is the $N$-point twiddle factor.

### 2.3 Special Cases: DC and Nyquist

At $k = 0$ and $k = N/2$, the formula simplifies. Both bins are purely real and can be recovered from the single complex value $Z[0]$:

$$X[0] = \text{Re}(Z[0]) + \text{Im}(Z[0])$$

$$X\!\left[\tfrac{N}{2}\right] = \text{Re}(Z[0]) - \text{Im}(Z[0])$$

This means DC and Nyquist share one physical slot in the packed buffer.

### 2.4 Hermitian Symmetry

For $k > N/2$, the spectrum is obtained by conjugation:

$$X[k] = \overline{X[N - k]}, \quad k = \frac{N}{2}+1, \ldots, N-1$$

Only $N/2 + 1$ unique spectral bins exist. The entire $N$-point spectrum can be reconstructed from bins $0$ through $N/2$.

### 2.5 Inverse: Repacking Spectral Bins Back to $Z[k]$

Given modified spectral bins $X[k]$ and $X[N/2 - k]$, we need to recover the packed complex values $Z[k]$ and $Z[N/2-k]$ so that a standard $N/2$-point IFFT can reconstruct the time-domain signal.

From the forward relation $X[k] = A_k - j B_k W_N^k$, we can extract:

$$A_k = \frac{X[k] + \overline{X[\tfrac{N}{2}-k]}}{2}$$

$$D_k = \frac{X[k] - \overline{X[\tfrac{N}{2}-k]}}{2} = -j B_k W_N^k$$

Solving for $B_k$:

$$B_k = j \cdot D_k \cdot \overline{W_N^k}$$

And finally recovering the packed buffer:

$$\boxed{Z[k] = A_k + B_k, \qquad Z\!\left[\tfrac{N}{2}-k\right] = \overline{A_k - B_k}}$$

For the DC/Nyquist slot:

$$Z[0] = \frac{X[0] + X[\tfrac{N}{2}]}{2} + j \cdot \frac{X[0] - X[\tfrac{N}{2}]}{2}$$

## 3. Implementation

### 3.1 Class Overview

The implementation consists of two classes:

| Class | Purpose |
|-------|---------|  
| `sfft::RFFT<T, CplxT, U, N>` | Executes both forward and inverse Real FFT |
| `sfft::RFFT_View<T, CplxT, U, N>` | Provides unpacked spectral access and in-place transform |

Template parameters:
- `T` — scalar real type (e.g., `double`, `Q23`)
- `CplxT` — complex type for the data (e.g., `std::complex<double>`, `Q23Complex`)
- `U` — complex type for twiddle factors (e.g., `std::complex<double>`, `Q31Complex`)
- `N` — FFT size (power of 2)

### 3.2 Forward Transform: `RFFT::process`

```cpp
static auto process(T *data)
{
     auto cdata = InterleavedComplexView<T, CplxT>(data);

    using TwidGen = TwiddleGenerator<U, N>;
    using StrideGen = StridedTwiddleGenerator<TwidGen, 2>;
    auto fft = FFT<CplxT, U, N / 2, StrideGen>();
    fft.process(cdata);

    return RFFT_View<T, CplxT, U, N>(data);
}
```

**Steps:**

1. **Pack**: Reinterpret the $N$-element real array as an $N/2$-element complex array (zero-cost).
2. **FFT**: Compute the $N/2$-point complex DIF FFT. The twiddle generator is strided by 2, which selects every other twiddle from the $N$-point table — exactly the twiddles needed for an $N/2$-point FFT.
3. **Return View**: Instead of materializing the full spectrum, return an `RFFT_View` that lazily unpacks bins on demand.

After this call, the buffer contains $Z[k]$ in **bit-reversed order** (the natural output of a DIF FFT). `RFFT_View` wraps the packed buffer in a `BitReversedView`, so `operator[]` uses natural indices regardless of physical layout.

### 3.3 Spectral Access: `RFFT_View::operator[]`

The view provides $O(1)$ access to any bin $X[k]$ for $k = 0, \ldots, N-1$:

```cpp
auto operator[](size_t k) const
```

**Logic:**

1. If $k \geq N/2$: set `wrap = true`, replace $k \leftarrow N - k$ (Hermitian symmetry).
2. If $k = 0$: return $(\text{Re}(Z[0]) + \text{Im}(Z[0]),\; 0)$.
3. If $k = N/2$: return $(\text{Re}(Z[0]) - \text{Im}(Z[0]),\; 0)$.
4. Otherwise, compute $A_k$, $B_k$, and apply the twiddle rotation:
   $$X[k] = A_k + (-j \cdot B_k) \cdot W_N^k$$
5. If `wrap`: return $\overline{X[k]}$.

The `BitReversedView` handles the bit-reversal indexing transparently — `_buffer[k]` returns $Z[k]$ regardless of the physical memory layout.

### 3.4 In-Place Spectral Manipulation: `RFFT_View::transform`

```cpp
template <typename Func>
void transform(Func &&f)
```

This method allows modifying the spectrum in-place without explicitly unpacking and repacking. It iterates over the unique bin pairs $(k, N/2-k)$ for $k = 0, \ldots, N/4$:

1. **Unpack** both bins $X[k]$ and $X[N/2-k]$ via `operator[]`.
2. **Apply** the user-provided function: `new_k = f(k, val_k)`.
3. **Repack** the modified values back into the packed buffer via `set_internal_packed`.

#### Example: Low-Pass Filter

```cpp
using namespace sfft;
const size_t N = 256;
using ComplexType = std::complex<double>;
using RealType = ComplexType::value_type;

auto rfft = RFFT<RealType, ComplexType, ComplexType, N>();
auto view = rfft.process(buffer);

view.transform([](size_t k, ComplexType val) {
    if (k > 5) return ComplexType(0.0, 0.0);
    return val;
});

rfft.inverse(buffer);
```

### 3.5 Repacking: `set_internal_packed`

This private method writes modified spectral bins back into the packed complex buffer. It inverts the algebra from Section 2.5:

```
Ak = (new_Xk + conj(new_Xnk)) / 2
Dk = (new_Xk - conj(new_Xnk)) / 2
Bk = j * Dk * conj(Wk)

Z[k]     = Ak + Bk
Z[N/2-k] = conj(Ak - Bk)
```

Both `_buffer[k]` and `_buffer[N/2-k]` are written simultaneously, since they are algebraically coupled.

For the DC/Nyquist case ($k = 0$):

```
Z[0] = ((DC + Nyquist) / 2) + j * ((DC - Nyquist) / 2)
```

### 3.6 Inverse Transform: `RFFT::inverse`

```cpp
static T *inverse(T *data)
{
     auto cdata = InterleavedComplexView<T, CplxT>(data);

    using TwidGen = TwiddleGenerator<U, N>;
    using StrideGen = StridedTwiddleGenerator<TwidGen, 2>;
    auto fft = FFT<CplxT, U, N / 2, StrideGen>();
    fft.inverse(cdata);

    return data;
}
```

The inverse assumes the buffer contains valid packed complex data $Z[k]$ (in bit-reversed order). It runs the $N/2$-point complex DIT IFFT using the `inverse()` method of the FFT class, which produces $z[n]$ in natural order. Since $z[n] = x[2n] + j \cdot x[2n+1]$, the real array is automatically reconstructed in-place. The `inverse()` path also applies the $1/N$ scaling via per-stage halving.

## 4. Data Flow Diagram

```
Real input x[0..N-1]
        │
        ▼
     InterleavedComplexView ──► z[0..N/2-1]  (complex, zero-cost)
        │
        ▼
   N/2-point DIF FFT ──► Z[k]  (bit-reversed)
        │
        ▼
   RFFT_View  ◄── lazy unpacking via operator[]
   │        │
   │  transform(f)  ──► modify spectrum in-place
   │        │
   │        ▼
   │   set_internal_packed  ──► repack Z[k]
   │
   ▼
   N/2-point DIT IFFT (via inverse()) ──► z[n]  (natural order)
        │
        ▼
   Real output x[0..N-1]  (in-place)
```

## 5. Complexity

| Operation | Time | Memory |
|-----------|------|--------|
| Forward RFFT | $O(\frac{N}{2} \log \frac{N}{2})$ | In-place |
| `operator[k]` | $O(1)$ | No allocation |
| `transform(f)` | $O(N)$ | In-place |
| Inverse RFFT | $O(\frac{N}{2} \log \frac{N}{2})$ | In-place |

The entire pipeline operates **in-place** on the original real buffer — no auxiliary memory is allocated.
