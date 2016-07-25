/**
 * << detailed description >>
 *
 * @file ParticleIndexPair.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */

#ifndef READDY_MAIN_PARTICLEINDEXPAIR_H
#define READDY_MAIN_PARTICLEINDEXPAIR_H

#include <cstddef>
#include <stdexcept>
#include <ostream>
#include <boost/functional/hash.hpp>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace model {
                struct ParticleIndexPair {
                    std::size_t idx1, idx2;
                    ParticleIndexPair(std::size_t idx1, std::size_t idx2){
                        if(idx1 < idx2) {
                            ParticleIndexPair::idx1 = idx1;
                            ParticleIndexPair::idx2 = idx2;
                        } else if(idx1 > idx2){
                            ParticleIndexPair::idx1 = idx2;
                            ParticleIndexPair::idx2 = idx1;
                        } else {
                            throw std::runtime_error("pair must not have equal indices");
                        }
                    }

                    friend std::size_t hash_value(const ParticleIndexPair &pip) {
                        std::size_t seed = 0;
                        boost::hash_combine(seed, pip.idx1);
                        boost::hash_combine(seed, pip.idx2);
                        return seed;
                    }

                    friend bool operator==(const ParticleIndexPair &pip1, const ParticleIndexPair &pip2) {
                        return pip1.idx1 == pip2.idx1 && pip1.idx2 == pip2.idx2;
                    }

                    friend std::ostream& operator<<(std::ostream& os, const ParticleIndexPair &pip) {
                        os << "ParticleIndexPair(" << pip.idx1 << ", " << pip.idx2 <<")";
                        return os;
                    }
                };

                struct ParticleIndexPairHasher {
                    std::size_t operator()(const ParticleIndexPair &pip) const {
                        return hash_value(pip);
                    }
                };
            }
        }
    }
}
#endif //READDY_MAIN_PARTICLEINDEXPAIR_H
