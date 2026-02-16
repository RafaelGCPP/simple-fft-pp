#pragma once
#include <cstddef>
#include <algorithm> // for std::swap
#include <bit>       // for std::countr_zero

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

    template <size_t N>
    struct BitReversalPermutation
    {
        static constexpr std::array<size_t, N> indices = []()
        {
            std::array<size_t, N> arr{};
            constexpr size_t bits = std::countr_zero(N);
            for (size_t i = 0; i < N; ++i)
            {
                arr[i] = reverse_bits<bits>(i);
            }
            return arr;
        }();
    };

    /**
     * @brief Decorator that provides a linear view and an in-place commit method.
     */
    template <typename T, size_t N, typename View = T *>
    class BitReversedView
    {
        View _data;
        static constexpr size_t _bits = std::countr_zero(N);

    public:
        constexpr explicit BitReversedView(View data) : _data(data) {}

        // Accessor (The "Decorator" part)
        constexpr decltype(auto) operator[](size_t index)
        {
            return _data[BitReversalPermutation<N>::indices[index]];
        }

        constexpr decltype(auto) operator[](size_t index) const
        {
            return _data[BitReversalPermutation<N>::indices[index]];
        }

        /**
         * @brief "Commits" the bit-reversal by physically reordering the data in RAM.
         * This makes the buffer linear so standard pointers/loops work.
         */
        constexpr void commit()
        {
            for (size_t i = 0; i < N; ++i)
            {
                size_t j = BitReversalPermutation<N>::indices[i];
                if (i < j)
                {
                    T temp = _data[i];
                    _data[i] = _data[j];
                    _data[j] = temp;
                }
            }
        }

        /**
         * @brief Apply a transform function to all elements in-place.
         * @param f Function taking (size_t index, T value) and returning T
         *
         * Example:
         * @code
         * view.transform([](size_t k, Complex val) {
         *     if (k > 20) return Complex(0.0, 0.0);
         *     return val;
         * });
         * @endcode
         */
        template <typename Func>
        constexpr void transform(Func &&f)
        {
            for (size_t i = 0; i < N; ++i)
            {
                (*this)[i] = f(i, (*this)[i]);
            }
        }

        static constexpr size_t size() { return N; }
    };

    template <typename T, typename CplxT>
    struct InterleavedComplexView
    {
        T *_data;
        constexpr explicit InterleavedComplexView(T *data) : _data(data) {}

        struct Proxy
        {
            T *_ptr;
            constexpr operator CplxT() const
            {
                return CplxT(_ptr[0], _ptr[1]);
            }

            constexpr Proxy &operator=(const CplxT &c)
            {
                _ptr[0] = c.real();
                _ptr[1] = c.imag();
                return *this;
            }

            constexpr Proxy &operator=(const Proxy &other)
            {
                if (this != &other)
                {
                    _ptr[0] = other._ptr[0];
                    _ptr[1] = other._ptr[1];
                }
                return *this;
            }

            constexpr T real() const { return _ptr[0]; }
            constexpr T imag() const { return _ptr[1]; }
        };

        constexpr Proxy operator[](size_t i) const
        {
            return Proxy{_data + 2 * i};
        }
    };

}