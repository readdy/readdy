#include <readdy/common/RandomProvider.h>
#include <boost/random.hpp>
#include <readdy/common/make_unique.h>

/**
 * << detailed description >>
 *
 * @file Random.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */
struct readdy::utils::RandomProvider::Impl {
    std::shared_ptr<boost::random::ranlux3> gen = std::make_shared<boost::random::ranlux3>();
};

double readdy::utils::RandomProvider::getNormal(double mean, double variance) {
    boost::random::normal_distribution<> n(mean, variance);
    return n(*pimpl->gen);
}

readdy::utils::RandomProvider::RandomProvider() : pimpl(std::make_unique<Impl>()) {

}

readdy::utils::RandomProvider &readdy::utils::RandomProvider::operator=(RandomProvider &&rhs) = default;

readdy::utils::RandomProvider::RandomProvider(RandomProvider &&rhs) = default;

readdy::utils::RandomProvider::~RandomProvider() = default;
