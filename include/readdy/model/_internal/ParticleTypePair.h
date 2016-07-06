/**
 * This header file contains the ParticleTypePair class. Its purpose is to provide a _sorted_ tuple of particle types,
 * so that it does not matter if, e.g., in a map, the type pair <i,j> or type pair <j,i> is requested.
 *
 * @file ParticleTypePair.h
 * @brief Header file containing the ParticleTypePair class, which is a sorted tuple of unsigned ints.
 * @author clonker
 * @date 09.06.16
 */

#ifndef READDY_MAIN_PARTICLETYPEPAIR_H
#define READDY_MAIN_PARTICLETYPEPAIR_H

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
