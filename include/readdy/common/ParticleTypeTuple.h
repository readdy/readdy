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
 * Definitions of hashers and equality operators of particle type pairs, triples, quadruples.
 *
 * @file ParticleTypeTuple.h
 * @brief Definitions of everything particle type tuple related
 * @author clonker
 * @date 20.03.17
 * @copyright BSD-3
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
            for_each_in_tuple_reverse(tuple, [&seed](ParticleTypeId ptt) { hash::combine(seed, ptt); });
        } else {
            for_each_in_tuple(tuple, [&seed](ParticleTypeId ptt) { hash::combine(seed, ptt); });
        }
        return seed;
    }
};

template<typename T>
class ForwardTupleHasher {
public:
    std::size_t operator()(const T& tuple) const {
        std::size_t seed {0};
        for_each_in_tuple(tuple, [&seed](ParticleTypeId ptt) { hash::combine(seed, ptt); });
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

using particle_type_pair = std::tuple<ParticleTypeId, ParticleTypeId>;
using particle_type_pair_hasher = ForwardBackwardTupleHasher<particle_type_pair>;
using particle_type_pair_equal_to = ForwardBackwardTupleEquality<particle_type_pair>;
template<typename T> using particle_type_pair_unordered_map = std::unordered_map<particle_type_pair, T, particle_type_pair_hasher, particle_type_pair_equal_to>;
using particle_type_triple = std::tuple<ParticleTypeId, ParticleTypeId, ParticleTypeId>;
using particle_type_triple_hasher = ForwardBackwardTupleHasher<particle_type_triple>;
using particle_type_triple_equal_to = ForwardBackwardTupleEquality<particle_type_triple>;
template<typename T> using particle_type_triple_unordered_map = std::unordered_map<particle_type_triple, T, particle_type_triple_hasher, particle_type_triple_equal_to>;
using particle_type_quadruple = std::tuple<ParticleTypeId, ParticleTypeId, ParticleTypeId, ParticleTypeId>;
using particle_type_quadruple_hasher = ForwardBackwardTupleHasher<particle_type_quadruple>;
using particle_type_quadruple_equal_to = ForwardBackwardTupleEquality<particle_type_quadruple>;
template<typename T> using particle_type_quadruple_unordered_map = std::unordered_map<particle_type_quadruple, T, particle_type_quadruple_hasher, particle_type_quadruple_equal_to>;

inline particle_type_triple sortTypeTriple(ParticleTypeId t1, ParticleTypeId t2, ParticleTypeId t3) {
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
sortTypeQuadruple(ParticleTypeId t1, ParticleTypeId t2, ParticleTypeId t3, ParticleTypeId t4) {
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
