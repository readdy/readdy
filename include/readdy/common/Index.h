/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * This header file contains definitions for the Index class, which makes the access to contiguously stored multi-
 * dimensional arrays easier.
 *
 * @file Index.h
 * @brief Definitions of the Index class, making contiguously stored multi-dimensional array access easier.
 * @author clonker
 * @date 11.09.17
 * @copyright BSD-3
 */

#pragma once

#include <initializer_list>
#include <array>
#include <numeric>
#include <tuple>
#include <cmath>
#include "tuple_utils.h"

namespace readdy::util {

namespace detail {
template<typename T1, typename... T>
struct variadic_first {
    /**
     * type of the first element of a variadic type tuple
     */
    using type = typename std::decay<T1>::type;
};
template<bool B, typename T = void> using disable_if = std::enable_if<!B, T>;

template<std::size_t N, typename T, typename Tuple>
struct all_of_type_impl {
    /**
     * type of the n-th element of a packed variadic type tuple
     */
    using nth_elem = typename std::decay<typename std::tuple_element<N - 1, Tuple>::type>::type;
    /**
     * if this is true, all elements in the parameter pack were of type T
     */
    static constexpr bool value = std::is_same<nth_elem, T>::value && all_of_type_impl<N - 1, T, Tuple>::value;
};

template<typename T, typename Tuple>
struct all_of_type_impl<0, T, Tuple> {
    /**
     * this is true for empty parameter packs
     */
    static constexpr bool value = true;
};

template<typename T, typename ...Args>
struct all_of_type {
    /**
     * if this is true, all elements of the parameter pack are of type T
     */
    static constexpr bool value = all_of_type_impl<sizeof...(Args), T, typename std::tuple<Args...>>::value;
};

}

namespace detail {
template<typename... Ix>
struct ComputeIndex {
    template<typename Strides, typename Indices = std::make_index_sequence<sizeof...(Ix)>>
    static constexpr auto compute(const Strides &strides, Ix &&... ix) {
        std::tuple<Ix...> tup(std::forward<Ix>(ix)...);
        return compute(strides, tup, Indices{});
    }

    template<typename Arr, std::size_t... I>
    static constexpr auto compute(const Arr &strides, const std::tuple<Ix...> &tup, std::index_sequence<I...>) {
        return (0 + ... + (strides[I] * std::get<I>(tup)));;
    }
};

}

template<std::size_t Dims, typename T = std::array<std::size_t, Dims>>
class Index {
    static_assert(Dims > 0, "Dims has to be > 0");
public:
    using GridDims = T;

    template<typename Container>
    static auto make_index(const Container& container) {
        return make_index(begin(container), end(container));
    }

    template<typename It>
    static auto make_index(It shapeBegin, It shapeEnd) {
        GridDims dims;
        std::copy(shapeBegin, shapeEnd, begin(dims));
        auto n_elems = std::accumulate(begin(dims), end(dims), 1u, std::multiplies<value_type>());

        GridDims strides;
        strides[0] = n_elems / dims[0];
        for (std::size_t d = 0; d < Dims - 1; ++d) {
            strides[d + 1] = strides[d] / dims[d + 1];
        }

        return Index<Dims, GridDims>{dims, strides, n_elems};
    }

    /**
     * The value type, inherited from GridDims::value_type
     */
    using value_type = typename GridDims::value_type;

    /**
     * Constructs an empty index object of specified dimensionality. Not of much use, really.
     */
    Index() : _size(), _cum_size(), n_elems(0) {}

    template<typename Shape>
    Index(const Shape &size)
            : _size(), n_elems(std::accumulate(begin(size), end(size), 1u, std::multiplies<value_type>())) {
        std::copy(begin(size), end(size), begin(_size));

        GridDims strides;
        strides[0] = n_elems / size[0];
        for (std::size_t d = 0; d < Dims - 1; ++d) {
            strides[d + 1] = strides[d] / size[d + 1];
        }
        _cum_size = std::move(strides);
    }

    /**
     * Constructs an index object with a number of size_t arguments that must coincide with the number of dimensions,
     * specifying the grid.
     * @tparam Args the argument types, must all be size_t
     * @param args the arguments
     */
    Index(GridDims size, GridDims strides, value_type nElems) : _size(std::move(size)), _cum_size(std::move(strides)),
                                                                n_elems(nElems) {}

    /**
     * the number of elements in this index, exactly the product of the grid dimensions
     * @return the number of elements
     */
    value_type size() const {
        return n_elems;
    }

    /**
     * the number of elements in this index
     * @return the number of elements
     */
    value_type nElements() const {
        return n_elems;
    }

    /**
     * Retrieve size of N-th axis
     * @tparam N the axis
     * @return size of N-th axis
     */
    template<int N>
    constexpr value_type get() const {
        return _size[N];
    }

    /**
     * retrieve size of N-th axis
     * @param N N
     * @return size of N-th axis
     */
    template<typename D>
    constexpr value_type operator[](D N) const {
        return _size[N];
    }

    /**
     * map Dims-dimensional index to 1D index
     * @tparam Ix the d-dimensional index template param type
     * @param ix the d-dimensional index
     * @return the 1D index
     */
    template<typename... Ix, typename Indices = std::make_index_sequence<Dims>>
    constexpr value_type operator()(Ix &&... ix) const {
        static_assert(std::size_t(sizeof...(ix)) == Dims, "wrong input dim");
        return detail::ComputeIndex<Ix...>::compute(_cum_size, std::forward<Ix>(ix)...);
    }

    /**
     * map Dims-dimensional array to 1D index
     * @param indices the Dims-dimensional index
     * @return the 1D index
     */
    template<typename Arr>
    value_type index(const Arr &indices) const {
        std::size_t result{0};
        for (std::size_t i = 0; i < Dims; ++i) {
            result += _cum_size[i] * indices[i];
        }
        return result;
    }

    /**
     * Inverse mapping 1D index to Dims-dimensional tuple
     * @param idx
     * @return
     */
    GridDims inverse(std::size_t idx) const {
        GridDims result;
        auto prefactor = n_elems / _size[0];
        for (std::size_t d = 0; d < Dims - 1; ++d) {
            auto x = std::floor(idx / prefactor);
            result[d] = x;
            idx -= x * prefactor;
            prefactor /= _size[d + 1];
        }
        result[Dims - 1] = idx;
        return result;
    }

private:
    GridDims _size;
    GridDims _cum_size;
    value_type n_elems;
};

/**
 * Definition of a 1D index object, basically an identity operator
 */
using Index1D = Index<1>;
/**
 * Definition of a 2D index object
 */
using Index2D = Index<2>;
/**
 * Definition of a 3D index object
 */
using Index3D = Index<3>;

}
