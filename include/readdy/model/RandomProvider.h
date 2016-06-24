/**
 * << detailed description >>
 *
 * @file Random.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_RANDOMPROVIDER_H
#define READDY_MAIN_RANDOMPROVIDER_H

#include <memory>
#include "Vec3.h"

namespace readdy {
    namespace model {
        class RandomProvider {
        public:
            RandomProvider();
            ~RandomProvider();
            RandomProvider(RandomProvider&& rhs);
            RandomProvider& operator=(RandomProvider&& rhs);
            virtual double getNormal(const double &mean = 0.0, const double &variance = 1.0);
            virtual double getUniform(double a = 0.0, double b = 1.0);
            virtual Vec3 getNormal3 (const double &mean = 0.0, const double &variance = 1.0);
        private:
            struct Impl;
            std::unique_ptr<Impl> pimpl;
        };
    }
}
#endif //READDY_MAIN_RANDOMPROVIDER_H
