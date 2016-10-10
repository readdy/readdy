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
    static short counter;

    const std::string name;
    const int order;
    const short id;

public:
    Potential(const std::string &name, const int order) : name(name), order(order), id(counter++) { }

    virtual ~Potential() = default;

    const short getId() const {
        return id;
    }

    const std::string &getName() const {
        return name;
    }

    const int getOrder() const {
        return order;
    }

    virtual double getMaximalForce(double kbt) const noexcept = 0;

};

}
}
}

#endif //READDY_MAIN_POTENTIAL_H
