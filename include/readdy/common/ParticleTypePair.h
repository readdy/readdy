/**
 * << detailed description >>
 *
 * @file ParticleTypePair.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */

#ifndef READDY_MAIN_PARTICLETYPEPAIR_H
#define READDY_MAIN_PARTICLETYPEPAIR_H

#include <cstddef>
#include <tuple>
#include <boost/functional/hash.hpp>

namespace readdy {
    namespace util {
        struct ParticleTypePair {
            unsigned int t1, t2;

            ParticleTypePair(unsigned int t1, unsigned int t2) {
                if (t1 <= t2) {
                    ParticleTypePair::t1 = t1;
                    ParticleTypePair::t2 = t2;
                } else {
                    ParticleTypePair::t1 = t2;
                    ParticleTypePair::t2 = t1;
                }
            }

            friend std::size_t hash_value(const ParticleTypePair &pair) {
                std::size_t seed = 0;
                boost::hash_combine(seed, pair.t1);
                boost::hash_combine(seed, pair.t2);
                return seed;
            }

            friend bool operator==(const ParticleTypePair &p1, const ParticleTypePair &p2) {
                return p1.t1 == p2.t1 && p1.t2 == p2.t2;
            }
        };

        class ParticleTypePairHasher {
        public:
            std::size_t operator()(const ParticleTypePair &k) const {
                return hash_value(k);
            }

            std::size_t operator()(const std::tuple<unsigned int, unsigned int> &k) const {
                std::size_t seed = 0;
                const auto &t1 = std::get<0>(k);
                const auto &t2 = std::get<1>(k);
                if (t1 <= t2) {
                    boost::hash_combine(seed, t1);
                    boost::hash_combine(seed, t2);
                } else {
                    boost::hash_combine(seed, t2);
                    boost::hash_combine(seed, t1);
                }
                return seed;
            }
        };
    }
}
#endif //READDY_MAIN_PARTICLETYPEPAIR_H
