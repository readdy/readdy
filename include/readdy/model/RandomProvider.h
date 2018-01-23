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
