#include <readdy/common/RandomProvider.h>
#include <boost/random.hpp>

/**
 * << detailed description >>
 *
 * @file Random.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 * @todo use pimpl idiom instead of static (mpi)
 */


double readdy::utils::RandomProvider::getNormal(double mean, double variance) {
    static boost::random::ranlux3 gen;
    boost::random::normal_distribution<> n (mean, variance);
    return n(gen);
}

