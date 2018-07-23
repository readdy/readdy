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
 * The random provider can provide normal and uniform distributed random numbers. The choice of random generator
 * can be altered by template parameter. Current default: mt19937.
 *
 * @file RandomProvider.h
 * @brief Header file containing the definitions for readdy::model::RandomProvider.
 * @author clonker
 * @date 19.04.16
 */

#pragma once
#include <memory>
#include <random>
#include <ctime>
#include <thread>
#include "readdy/common/common.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(rnd)

template<typename RealType=scalar, typename Generator = std::default_random_engine>
RealType normal(const RealType mean = 0.0, const RealType variance = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::normal_distribution<RealType> distribution(mean, variance);
    return distribution(generator);
}

template<typename RealType=scalar, typename Generator = std::default_random_engine>
RealType uniform_real(const RealType a = 0.0, const RealType b = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_real_distribution<RealType> distribution(a, b);
    return distribution(generator);
}

template<typename IntType=int, typename Generator = std::default_random_engine>
IntType uniform_int(const IntType a, const IntType b) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_int_distribution<IntType> distribution(a, b);
    return distribution(generator);
}

template<typename RealType=scalar, typename Generator = std::default_random_engine>
RealType exponential(RealType lambda = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::exponential_distribution<RealType> distribution(lambda);
    return distribution(generator);
}

template<typename scalar, typename Generator = std::default_random_engine>
Vec3 normal3(const scalar mean = 0.0, const scalar variance = 1.0) {
    return {normal<scalar, Generator>(mean, variance),
            normal<scalar, Generator>(mean, variance),
            normal<scalar, Generator>(mean, variance)};
}

template<typename Iter, typename Gen = std::default_random_engine>
Iter random_element(Iter start, const Iter end) {
    using IntType = typename std::iterator_traits<Iter>::difference_type;
    std::advance(start, uniform_int<IntType, Gen>(0, std::distance(start, end)));
    return start;
}

NAMESPACE_END(rnd)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
