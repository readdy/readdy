/**
 * << detailed description >>
 *
 * @file ParticleTypePair.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#ifndef READDY_MAIN_PARTICLETYPEPAIR_H
#define READDY_MAIN_PARTICLETYPEPAIR_H

#include <boost/mpl/size_t.hpp>
#include <boost/functional/hash.hpp>

namespace readdy {
    namespace model {
        namespace _internal {
            struct ParticleTypePair {
                unsigned int t1, t2;

                ParticleTypePair(unsigned int t1, unsigned int t2) {
                    if(t1 <= t2) {
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
        }
    }
}
#endif //READDY_MAIN_PARTICLETYPEPAIR_H
