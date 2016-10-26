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

template<typename Generator = std::default_random_engine>
double normal(const double mean = 0.0, const double variance = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::normal_distribution<double> distribution(mean, variance);
    return distribution(generator);
}

template<typename Generator = std::default_random_engine>
double uniform_real(const double a = 0.0, const double b = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_real_distribution<double> distribution(a,b);
    return distribution(generator);
}

template<typename Generator = std::default_random_engine>
double uniform_int(const double a = 0.0, const double b = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_int_distribution<double> distribution(a,b);
    return distribution(generator);
}

template<typename Generator = std::default_random_engine>
double exponential(double lambda = 1.0) {
    static thread_local Generator generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::exponential_distribution<double> distribution(lambda);
    return distribution(generator);
}

template<typename Generator = std::default_random_engine>
Vec3 normal3(const double mean = 0.0, const double variance = 1.0) {
    return {normal<Generator>(mean, variance), normal<Generator>(mean, variance), normal<Generator>(mean, variance)};
}

template <typename Iter, typename Generator = std::default_random_engine>
Iter random_element(Iter start, Iter end) {
    std::advance(start, uniform_int<Generator>(0, std::distance(start, end)));
    return start;
}


}
}
}
#endif //READDY_MAIN_RANDOMPROVIDER_H
