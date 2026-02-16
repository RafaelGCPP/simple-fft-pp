#pragma once
#include "view_types.h"
#include "twiddles.h"

namespace sfft
{

    template <typename T>
    constexpr T scale_in_half(const T &value)
    {
        return value * T(0.5);
    }

    template <typename T, typename U, int N, typename TwidGen = TwiddleGenerator<U, N>>
    class FFT
    {

    public:
        FFT() = default;

        template <typename View>
        constexpr auto process(View data)
        {
            int stride, num_blocks;
            for (stride = N / 2, num_blocks = 1; stride >= 1; stride /= 2, num_blocks *= 2)
            {
                for (int block = 0; block < num_blocks; ++block)
                {
                    int block_start = block * 2 * stride;
                    for (int j = 0; j < stride; ++j)
                    {
                        T a = data[block_start + j];
                        T b = data[block_start + j + stride];
                        U twiddle = TwidGen::twiddles[j * num_blocks];

                        data[block_start + j] = a + b;
                        T diff = a - b;
                        data[block_start + j + stride] = diff * twiddle;
                    }
                }
            }
            return BitReversedView<T, N, View>(data); // Return a view for bit-reversed access
        }

        template <typename View>
        constexpr auto inverse(View data)
        {
            int stride, num_blocks;
            for (stride = 1, num_blocks = N / 2; stride < N; stride *= 2, num_blocks /= 2)
            {
                for (int block = 0; block < num_blocks; ++block)
                {
                    int block_start = block * 2 * stride;
                    for (int j = 0; j < stride; ++j)
                    {
                        T a = data[block_start + j];
                        T b = data[block_start + j + stride];
                        U twiddle = conj(TwidGen::twiddles[j * num_blocks]);
                        b *= twiddle;
                        data[block_start + j] = scale_in_half(a + b);
                        data[block_start + j + stride] = scale_in_half(a - b);
                    }
                }
            }
            return data;
        }
    };

} // namespace sfft
