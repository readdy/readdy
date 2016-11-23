/**
 * The random provider can provide normal and uniform distributed random numbers. The choice of random generator
 * can be altered by template parameter. Current default: mt19937.
 *
 * @file Random.h
 * @brief Header file containing the definitions for readdy::model::RandomProvider.
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_RANDOMPROVIDER_H
#define READDY_MAIN_RANDOMPROVIDER_H

#include <memory>
#include <random>
#include <time.h>
#include <thread>
#include "Vec3.h"

namespace readdy {
namespace model {
namespace rnd {

template<typename RealType=double, typename Generator = std::default_random_engine>
RealType normal(const RealType mean = 0.0, const RealType variance = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::normal_distribution<RealType> distribution(mean, variance);
    return distribution(generator);
}

template<typename RealType=double, typename Generator = std::default_random_engine>
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

template<typename RealType=double, typename Generator = std::default_random_engine>
RealType exponential(RealType lambda = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::exponential_distribution<RealType> distribution(lambda);
    return distribution(generator);
}

template<typename Generator = std::default_random_engine>
Vec3 normal3(const double mean = 0.0, const double variance = 1.0) {
    return {normal<double, Generator>(mean, variance),
            normal<double, Generator>(mean, variance),
            normal<double, Generator>(mean, variance)};
}

template<typename Iter, typename Gen = std::default_random_engine>
Iter random_element(Iter start, const Iter end) {
    using IntType = typename std::iterator_traits<Iter>::difference_type;
    std::advance(start, uniform_int<IntType, Gen>(0, std::distance(start, end)));
    return start;
}


}
}
}
#endif //READDY_MAIN_RANDOMPROVIDER_H
