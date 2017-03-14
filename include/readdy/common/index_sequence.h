/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
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
 * @file index_sequence.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.08.16
 */

#pragma once

#ifdef READDY_CPP14

#include <utility>

#else

#include <cstddef>

namespace std {
    template<size_t... Ints>
    struct index_sequence {
        using type = index_sequence;
        using value_type = size_t;

        static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
    };

    // --------------------------------------------------------------

    template<class Sequence1, class Sequence2>
    struct _merge_and_renumber;

    template<size_t... I1, size_t... I2>
    struct _merge_and_renumber<index_sequence<I1...>, index_sequence<I2...>>
            : index_sequence<I1..., (sizeof...(I1) + I2)...> {
    };

    // --------------------------------------------------------------

    template<size_t N>
    struct make_index_sequence
            : _merge_and_renumber<typename make_index_sequence<N / 2>::type,
                    typename make_index_sequence<N - N / 2>::type> {
    };

    template<>
    struct make_index_sequence<0> : index_sequence<> { };
    template<>
    struct make_index_sequence<1> : index_sequence<0> { };
}
#endif
