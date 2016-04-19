/**
 * << detailed description >>
 *
 * @file Random.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY2_MAIN_RANDOMPROVIDER_H
#define READDY2_MAIN_RANDOMPROVIDER_H

namespace readdy {
    namespace utils {
        class RandomProvider {
        public:
            virtual double getNormal(double mean = 0.0, double variance = 1.0);
        };
    }
}
#endif //READDY2_MAIN_RANDOMPROVIDER_H
