#pragma once

#include "fft_core.h"

template <typename T, typename CplxT, typename TwidGen>
class RFFT_View
{
    // Usamos o seu Decorator internamente
    static constexpr size_t N = TwidGen::N_Value;
    BitReversedView<CplxT, N/2> _buffer;

public:
    RFFT_View(T *data) : _buffer(reinterpret_cast<CplxT *>(data)) {}

    auto operator[](size_t k) const
    {
        // Agora o operator[] da BitReversedView faz o reverse_bits pra você!
        // Não importa se o buffer está bagunçado ou se você já deu 'commit()',
        // a lógica de descompactação (Ak, Bk) permanece linear.
        bool  wrap= (k >= N/2); // Se k está na segunda metade do espectro, precisamos acessar a simetria
        if (wrap) {
            k=N-k; // Explora a simetria para acessar a segunda metade do espectro
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

        return (wrap ? conj(result) : result); // Se k estava na segunda metade, retorna a simetria conjugada>;
    }
};

template <typename T, typename CplxT, typename TwidGen>
class RFFT
{
public:
    static auto process(T *data)
    {
        // constexpr size_t N = TwidGen::N_Value;

        // First reshape real input into packed complex.
        auto cdata = reinterpret_cast<CplxT *>(data);

        using StrideGen = StridedTwiddleGenerator<TwidGen, 2>;
        auto fft = FFT<CplxT, StrideGen>();
        fft.process(cdata);

        return RFFT_View<T, CplxT, TwidGen>(data); // Return a view for unpacked access
    }
};
