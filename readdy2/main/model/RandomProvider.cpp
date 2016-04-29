#include <readdy/model/RandomProvider.h>
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
namespace readdy {
    namespace model {
        struct RandomProvider::Impl {
            std::shared_ptr<boost::random::ranlux3> gen = std::make_shared<boost::random::ranlux3>();
        };

        double RandomProvider::getNormal(double mean, double variance) {
            boost::random::normal_distribution<> n(mean, variance);
            return n(*pimpl->gen);
        }

        RandomProvider::RandomProvider() : pimpl(std::make_unique<Impl>()) {

        }

        Vec3 RandomProvider::getNormal3(double mean, double variance) {
            return {getNormal(mean, variance), getNormal(mean, variance), getNormal(mean, variance)};
        }


        RandomProvider &RandomProvider::operator=(RandomProvider &&rhs) = default;

        RandomProvider::RandomProvider(RandomProvider &&rhs) = default;

        RandomProvider::~RandomProvider() = default;

    }
}
