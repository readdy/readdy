/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file Index.h
 * @brief << brief description >>
 * @author clonker
 * @date 11.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <initializer_list>
#include <array>
#include <numeric>
#include "macros.h"
#include "tuple_utils.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

template<std::size_t Dims>
class Index {
    static_assert(Dims > 0, "Dims has to be > 0");
public:
    using GridDims = std::array<int, Dims>;
    using value_type = typename GridDims::value_type;

    Index() : size(), n_elems(0) {}

    template<typename ...Args>
    Index(Args &&...args) : size({std::forward<Args>(args)...}),
                            n_elems(std::accumulate(size.begin(), size.end(), 1, std::multiplies<value_type>())) {}

    Index(GridDims &&arr) : size(std::move(arr)),
                            n_elems(std::accumulate(size.begin(), size.end(), 1, std::multiplies<value_type>())) {}

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
        return size[N];
    }

    /**
     * retrieve size of N-th axis
     * @param N N
     * @return size of N-th axis
     */
    constexpr value_type operator[](std::size_t N) const {
        return size[N];
    }

    /**
     * map Dims-dimensional index to 1D index
     * @tparam Ix the d-dimensional index template param type
     * @param index the d-dimensional index
     * @return the 1D index
     */
    template<typename... Ix>
    constexpr std::size_t operator()(Ix &&... index) const {
        static_assert(sizeof...(index) == Dims, "wrong input dim");
        return (*this)({std::forward<Ix>(index)...});

    }

    /**
     * map Dims-dimensional index to 1D index
     * @param ix the d-dimensional index
     * @return the 1D index
     */
    constexpr std::size_t operator()(const GridDims &ix) const {
        std::size_t result = 0;
        auto prefactor = n_elems / size[0];
        for(std::size_t d = 0; d < Dims-1; ++d) {
            result += prefactor * ix[d];
            prefactor /= size[d+1];
        }
        result += ix[Dims-1];
        return result;
    }

    /**
     * Inverse mapping 1D index to Dims-dimensional tuple
     * @param idx
     * @return
     */
    constexpr GridDims inverse(std::size_t idx) const {
        GridDims result;
        auto prefactor = n_elems / size[0];
        for(std::size_t d = 0; d < Dims-1; ++d) {
            auto x = std::floor(idx / prefactor);
            result[d] = x;
            idx -= x * prefactor;
            prefactor /= size[d+1];
        }
        result[Dims-1] = idx;
        return result;
    }

private:
    GridDims size;
    value_type n_elems;
};

using Index1D = Index<1>;
using Index2D = Index<2>;
using Index3D = Index<3>;

NAMESPACE_END(util)
NAMESPACE_END(readdy)