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
 * Definitions of hashers and equality operators of particle type pairs, triples, quadruples.
 *
 * @file ParticleTypeTuple.h
 * @brief Definitions of everything particle type tuple related
 * @author clonker
 * @date 20.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <utility>
#include <tuple>

#include <spdlog/fmt/ostr.h>

#include "common.h"
#include "hash.h"
#include "tuple_utils.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

template<typename T, std::size_t N = std::tuple_size<typename std::remove_reference<T>::type>::value>
class ForwardBackwardTupleHasher {
public:
    /**
     * Evaluates the hash of tuples independent of their reversedness.
     * @param tuple the tuple to compare
     * @return the hash of the reversed tuple if tuple[0] > tuple[-1], otherwise hash of the non-reversed version
     */
    std::size_t operator()(const T &tuple) const {
        std::size_t seed{0};
        // auto sorted_quadruple = call(sortTypeQuadruple, quadruple);
        if (std::get<0>(tuple) > std::get<N - 1>(tuple)) {
            for_each_in_tuple_reverse(tuple, [&seed](particle_type_type ptt) { hash::combine(seed, ptt); });
        } else {
            for_each_in_tuple(tuple, [&seed](particle_type_type ptt) { hash::combine(seed, ptt); });
        }
        return seed;
    }
};

template<typename T>
class ForwardTupleHasher {
public:
    std::size_t operator()(const T& tuple) const {
        std::size_t seed {0};
        for_each_in_tuple(tuple, [&seed](particle_type_type ptt) { hash::combine(seed, ptt); });
        return seed;
    }
};

template<typename T>
class ForwardTupleEquality {
public:
    bool operator()(const T &lhs, const T &rhs) const {
        return lhs == rhs;
    }
};


template<typename T>
class ForwardBackwardTupleEquality {
public:
    /**
     * Evaluates the equality of two tuples independently of their reversedness.
     * @param lhs the one tuple
     * @param rhs the other tuple
     * @return true if lhs == rhs or lhs == reverse(rhs)
     */
    bool operator()(const T &lhs, const T &rhs) const {
        return lhs == rhs || lhs == reverse(rhs);
    }
};

using particle_type_pair = std::tuple<particle_type_type, particle_type_type>;
using particle_type_pair_hasher = ForwardBackwardTupleHasher<particle_type_pair>;
using particle_type_pair_equal_to = ForwardBackwardTupleEquality<particle_type_pair>;
template<typename T> using particle_type_pair_unordered_map = std::unordered_map<particle_type_pair, T, particle_type_pair_hasher, particle_type_pair_equal_to>;
using particle_type_triple = std::tuple<particle_type_type, particle_type_type, particle_type_type>;
using particle_type_triple_hasher = ForwardBackwardTupleHasher<particle_type_triple>;
using particle_type_triple_equal_to = ForwardBackwardTupleEquality<particle_type_triple>;
template<typename T> using particle_type_triple_unordered_map = std::unordered_map<particle_type_triple, T, particle_type_triple_hasher, particle_type_triple_equal_to>;
using particle_type_quadruple = std::tuple<particle_type_type, particle_type_type, particle_type_type, particle_type_type>;
using particle_type_quadruple_hasher = ForwardBackwardTupleHasher<particle_type_quadruple>;
using particle_type_quadruple_equal_to = ForwardBackwardTupleEquality<particle_type_quadruple>;
template<typename T> using particle_type_quadruple_unordered_map = std::unordered_map<particle_type_quadruple, T, particle_type_quadruple_hasher, particle_type_quadruple_equal_to>;

inline particle_type_triple sortTypeTriple(particle_type_type t1, particle_type_type t2, particle_type_type t3) {
    if (t1 > t2) {
        std::swap(t1, t2);
    }
    if (t1 > t3) {
        std::swap(t1, t3);
    }
    if (t2 > t3) {
        std::swap(t2, t3);
    }
    return std::make_tuple(t1, t2, t3);
}

inline particle_type_quadruple
sortTypeQuadruple(particle_type_type t1, particle_type_type t2, particle_type_type t3, particle_type_type t4) {
    // divide into two sets {a, b} and {c, d} and sort them separately
    if (t1 > t2) {
        std::swap(t1, t2);
    }
    if (t3 > t4) {
        std::swap(t3, t4);
    }

    if (t1 > t3) {
        // t3 is globally smallest, compare (t1, t4)
        if (t1 > t4) {
            // t4 is next smallest, thus (t3, t4, t1, t2)
            return std::make_tuple(t3, t4, t1, t2);
        } else {
            // t1 is next smallest, thus (t3, t1, ...) -> compare (t2, t4)
            if (t2 > t4) {
                return std::make_tuple(t3, t1, t4, t2);
            } else {
                return std::make_tuple(t3, t1, t2, t4);
            }
        }
    } else {
        // t1 is globally smallest, compare (t2, t3)
        if (t2 > t3) {
            // t3 is next smallest, thus (t1, t3, ...) -> compare (t2, t4)
            if (t2 > t4) {
                return std::make_tuple(t1, t3, t4, t2);
            } else {
                return std::make_tuple(t1, t3, t2, t4);
            }
        } else {
            // t2 is next smallest, thus (t1, t2, t3, t4)
            return std::make_tuple(t1, t2, t3, t4);
        }
    }
}

NAMESPACE_END(util)
NAMESPACE_END(readdy)
