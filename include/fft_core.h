#pragma once
#include <cstddef>
#include <algorithm> // for std::swap

namespace sfft
{
    /**
     * @brief Compile-time bit reversal of an index
     */
    template <size_t Bits>
    constexpr size_t reverse_bits(size_t n)
    {
        size_t result = 0;
        for (size_t i = 0; i < Bits; ++i)
        {
            result = (result << 1) | (n & 1);
            n >>= 1;
        }
        return result;
    }

    /**
     * @brief Decorator that provides a linear view and an in-place commit method.
     */
    template <typename T, size_t N>
    class BitReversedView
    {
        T *const _data;
        static constexpr size_t _bits = std::countr_zero(N);

    public:
        explicit BitReversedView(T *data) : _data(data) {}

        // Accessor (The "Decorator" part)
        constexpr T &operator[](size_t index)
        {
            return _data[reverse_bits<_bits>(index)];
        }

        constexpr const T &operator[](size_t index) const
        {
            return _data[reverse_bits<_bits>(index)];
        }

        /**
         * @brief "Commits" the bit-reversal by physically reordering the data in RAM.
         * This makes the buffer linear so standard pointers/loops work.
         */
        void commit()
        {
            for (size_t i = 0; i < N; ++i)
            {
                size_t j = reverse_bits<_bits>(i);
                if (i < j)
                {
                    std::swap(_data[i], _data[j]);
                }
            }
        }

        static constexpr size_t size() { return N; }
    };

    template <typename T>
    constexpr T scale_in_half(const T &value)
    {
        return value * T(0.5);
    }

    template <typename T, typename U, int N, bool real_input = false>
    class FFT
    {
    private:
        U twiddles[N];

    public:
        FFT()
        {
            for (size_t i = 0; i < N / 2; ++i)
            {
                twiddles[i] = U(cos(-2.0 * M_PI * i / N), sin(-2.0 * M_PI * i / N));
            }
        }

        void process(T *data)
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
                        U twiddle = twiddles[j * num_blocks  ];

                        data[block_start + j] = a + b;
                        T diff = a - b;
                        data[block_start + j + stride] = diff * twiddle;
                    }
                }
            }
        }

        void inverse(T *data)
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
                        U twiddle = conj(twiddles[j * num_blocks]);
                        b = b * twiddle;
                        data[block_start + j] = scale_in_half(a + b);
                        data[block_start + j + stride] = scale_in_half(a - b);
                    }
                }
            }
        }
    };

} // namespace sfft
