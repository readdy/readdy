/**
 * This header file contains the base class for all potentials. In particular, all potentials shall have
 * a uuid by which they can be identified, a name that describes the type and an order, which indicates the number of
 * particles which are interacting with one another.
 *
 * @file Potential.h
 * @brief Header file containing the base class for all potentials.
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIAL_H
#define READDY_MAIN_POTENTIAL_H

#include <string>
#include <array>
#include <readdy/model/Particle.h>
#include <boost/uuid/random_generator.hpp>

namespace readdy {
    namespace model {
        namespace potentials {

            class Potential {
                const std::string name;
                const int order;
                boost::uuids::uuid id;

            public:
                Potential(const std::string &name, const int order) : name(name), order(order),
                                                                      id(boost::uuids::random_generator()()) {
                }
                virtual ~Potential() = default;

                const std::string &getName() const {
                    return name;
                }

                const int getOrder() const {
                    return order;
                }


                const boost::uuids::uuid &getId() const {
                    return id;
                }

                void setId(const boost::uuids::uuid id) {
                    Potential::id = id;
                }

                virtual Potential* replicate() const = 0;
            };

            namespace _internal {
                template<typename T>
                struct PotentialName { static const std::string value; };
            }

            template<typename PotentialType>
            const std::string& getPotentialName() {
                return readdy::model::potentials::_internal::PotentialName<PotentialType>::value;
            }
        }
    }
}

#endif //READDY_MAIN_POTENTIAL_H
