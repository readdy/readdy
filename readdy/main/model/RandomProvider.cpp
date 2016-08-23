/**
 * << detailed description >>
 *
 * @file Random.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include <readdy/model/RandomProvider.h>
#include <boost/random.hpp>
#include <readdy/common/make_unique.h>
#include <ctime>

namespace readdy {
    namespace model {
        struct RandomProvider::Impl {
            std::unique_ptr<boost::random::ranlux3> gen = std::make_unique<boost::random::ranlux3>((unsigned int) std::time(0));

            std::unique_ptr<boost::random::normal_distribution<>> normal01 = std::make_unique<boost::random::normal_distribution<>>(0.,1.);
            std::unique_ptr<boost::random::uniform_01<>> uniform01 = std::make_unique<boost::random::uniform_01<>>();
        };

        double RandomProvider::getNormal(const double mean, const double variance) {
            if(mean == 0 && variance == 1) return (*pimpl->normal01)(*pimpl->gen);
            boost::random::normal_distribution<> n(mean, variance);
            return n(*pimpl->gen);
        }

        RandomProvider::RandomProvider() : pimpl(std::make_unique<Impl>()) {
        }

        Vec3 RandomProvider::getNormal3(const double mean, const double variance) {
            return {getNormal(mean, variance), getNormal(mean, variance), getNormal(mean, variance)};
        }

        double RandomProvider::getUniform(double a, double b) {
            if(a > b) std::swap(a, b);
            return (*pimpl->uniform01)(*pimpl->gen)*(b-a) + a;
        }


        RandomProvider &RandomProvider::operator=(RandomProvider &&rhs) = default;

        RandomProvider::RandomProvider(RandomProvider &&rhs) = default;

        RandomProvider::~RandomProvider() = default;

    }
}
