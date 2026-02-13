#pragma once

#include "fft_core.h"

template <typename T, typename CplxT, typename TwidGen>
class RFFT_View
{
    // We use the Decorator internally
    static constexpr size_t N = TwidGen::N_Value;
    BitReversedView<CplxT, N / 2> _buffer;

public:
    RFFT_View(T *data) : _buffer(reinterpret_cast<CplxT *>(data)) {}

    auto operator[](size_t k) const
    {
        // The BitReversedView operator[] handles reverse_bits for you!
        // It doesn't matter if the buffer is scrambled or if you've already called 'commit()',
        // the unpacking logic (Ak, Bk) remains linear.
        bool wrap = (k >= N / 2); // If k is in the second half of the spectrum, we need to access the symmetry
        if (wrap)
        {
            k = N - k; // Exploit symmetry to access the second half of the spectrum
        }

        if (k == 0)
            return CplxT(_buffer[0].real() + _buffer[0].imag(), T(0));
        if (k == N / 2)
            return CplxT(_buffer[0].real() - _buffer[0].imag(), T(0));

        size_t knp = N / 2 - k;

        auto Ak = scale_in_half(_buffer[k] + conj(_buffer[knp]));
        auto Bk = scale_in_half(_buffer[k] - conj(_buffer[knp]));

        auto minus_j_Bk = CplxT(Bk.imag(), -Bk.real());
        auto W = TwidGen::get_twiddle(k);

        auto result = Ak + (minus_j_Bk * W);

        return (wrap ? conj(result) : result); // If k was in the second half, return the conjugate symmetry
    }
    template <typename Func>
    void transform(Func &&f)
    {
        // We only process the real half.
        // Each pair (k, N-k) modifies ONE slot in the physical buffer.
        for (size_t k = 0; k <= N / 4; ++k)
        {
            // 1. Get the current values (Unpack)
            CplxT val_k = (*this)[k];
            CplxT val_nk = (*this)[N/2 - k];

            // 2. Let the user apply their function
            CplxT new_k = f(k, val_k);
            CplxT new_nk = f(N/2 - k, val_nk);

            // 3. Repack (This is where the magic happens)
            // You need a set_packed(k, new_k, new_nk) function
            // that handles BitReverse + Symmetry Algebra.
            this->set_internal_packed(k, new_k, new_nk);
        }
    }

private:
    void set_internal_packed(size_t k, CplxT new_Xk, CplxT new_Xnk)
    {
        // Special case: DC and Nyquist share the same physical slot [0]
        if (k == 0)
        {
            T dc = new_Xk.real();
            T nyq = new_Xnk.real();
            _buffer[0] = CplxT(scale_in_half(dc + nyq), scale_in_half(dc - nyq));
            return;
        }

        size_t knp = N / 2 - k;

        // Recover the symmetric/anti-symmetric parts from the spectral bins:
        // Ak = (Xk + conj(X_{N/2-k})) / 2
        // Dk = (Xk - conj(X_{N/2-k})) / 2   (this is -j*Bk*Wk, NOT Bk itself)
        auto Ak = scale_in_half(new_Xk + conj(new_Xnk));
        auto Dk = scale_in_half(new_Xk - conj(new_Xnk));

        // To recover Bk from Dk:  Dk = -j*Bk*Wk  =>  Bk = j*Dk*conj(Wk)
        auto W = conj(TwidGen::get_twiddle(k));

        // Dk * conj(Wk)
        auto W_Dk = Dk * W;

        // j * W_Dk -> { -imag, real }
        auto j_W_Dk = CplxT(-W_Dk.imag(), W_Dk.real());

        // Bk = j * Dk * conj(Wk)
        // Z[k]     = Ak + Bk
        // Z[N/2-k] = conj(Ak - Bk)
        _buffer[k]   = Ak + j_W_Dk;
        _buffer[knp] = conj(Ak - j_W_Dk);
    }
};

template <typename T, typename CplxT, typename TwidGen>
class RFFT
{
public:
    static auto process(T *data)
    {
        // First reshape real input into packed complex.
        auto cdata = reinterpret_cast<CplxT *>(data);

        using StrideGen = StridedTwiddleGenerator<TwidGen, 2>;
        auto fft = FFT<CplxT, StrideGen>();
        fft.process(cdata);

        return RFFT_View<T, CplxT, TwidGen>(data); // Return a view for unpacked access
    }
};

template <typename T, typename CplxT, typename TwidGen>
class IRFFT
{
public:
    static T *process(T *data)
    {
        // First reshape real input into packed complex.
        auto cdata = reinterpret_cast<CplxT *>(data);

        using StrideGen = StridedTwiddleGenerator<TwidGen, 2>;
        auto ifft = IFFT<CplxT, StrideGen>();
        ifft.process(cdata);

        return data; // Return raw pointer for linear access after IRFFT
    }
};
