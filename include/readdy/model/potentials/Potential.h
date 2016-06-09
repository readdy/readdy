/**
 * << detailed description >>
 *
 * @file Potential.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIAL_H
#define READDY_MAIN_POTENTIAL_H

#include <string>
#include <array>
#include <readdy/model/Particle.h>

namespace readdy {
    namespace model {
        namespace potentials {

            class Potential {
                const std::string name;
                const int order;

            public:
                Potential(const std::string &name, const int order) : name(name), order(order) {
                }

                const std::string &getName() const {
                    return name;
                }

                const int getOrder() const {
                    return order;
                }

            };

            namespace _internal {
                template<typename T>
                struct PotentialName { static const std::string value; };

                template<typename PotentialType>
                const std::string& getPotentialName() {
                    return readdy::model::potentials::_internal::PotentialName<PotentialType>::value;
                }

            }
        }
    }
}
#endif //READDY_MAIN_POTENTIAL_H
